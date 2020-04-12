#include "commons.h"
namespace Search{
string Commons::getTimeAndDate() {
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(buffer,sizeof(buffer),"%d-%m-%Y-%Hh%Mm%Ss",timeinfo);
  std::string str(buffer);
  return str;
}

// I -> L, C+57.021 -> C(C), M+15.995 -> M(O)
string Commons::AdjustPeptide(const string& _peptide) {
  string peptide = _peptide;
  peptide = Utility::Commons::ReplaceIWithL(peptide);
  peptide = Utility::Commons::ReplacePTMOnCWithText(peptide, "(C)");
  peptide = Utility::Commons::ReplacePTMOnMWithText(peptide, "(O)");
  return peptide;
}

void Commons::ReadSpectraHelper(
    const HyperParams& params,
    const string& file, 
    vector<Spectrum*>* spectra, 
    unordered_map<int, vector<int>>* map_spectra_by_charge,
    unordered_map<int, vector<int>>* map_spectra_by_mass,
    unordered_map<string, vector<int>>* map_peptide_to_indices,
    unordered_map<string, int>* map_ms_title_to_index,
    bool peptide_i2l, bool peptide_ptm_replace) {
  int spectra_size = 0;
  bool verbose = false;
  Utility::IO::ReadSpectraFromMGF(spectra, 
                                  map_spectra_by_charge, 
                                  map_spectra_by_mass,
                                  map_peptide_to_indices,
                                  map_ms_title_to_index, 
                                  &spectra_size, 
                                  file, 
                                  params.mz_scale, 
                                  params.min_mz, 
                                  params.max_mz, 
                                  params.precision, 
                                  params.select_topk,
                                  params.window_mz, 
                                  static_cast<Utility::PEAK_INTENSITY_RESCALE_METHOD>(params.peak_intensity_rescale_method),
                                  params.remove_precursor,
                                  peptide_i2l,
                                  peptide_ptm_replace,
                                  params.filter_unfragmented_ms2);

  cout << "#spectra in total so far: " << (*spectra).size() << endl;
  
  cout << "STATS about the MGF files so far in total." << endl;
  for (const auto& entry : *map_spectra_by_charge) {
    cout << entry.first << " charge #spectra: " << entry.second.size() << endl; 
  }
  cout << endl;
}

void Commons::ReleaseMemoryHelper(
    string spectra_name, 
    vector<Spectrum*>* spectra) {
  Utility::Timer timer;
  cout << "Releasing memory of " << spectra_name << endl;
  for (int i = 0; i < (*spectra).size(); ++i) {
    delete (*spectra)[i];
  }
  cout << "Releasing memory takes: " << timer.stop() << endl;
  cout << endl;
}

void Commons::ReadMsgfIds(
    string fileName, 
    unordered_map<string, string>* map_title_to_peptide, 
    unordered_map<string, vector<string>>* map_peptide_to_titles) {

  map_title_to_peptide->clear();
  map_peptide_to_titles->clear();
  ifstream reader(fileName);
  if (!reader.is_open()) {
    cout << "Error! Can not open" << fileName << endl;
    exit(-1);
  }
  string line;
  // Get rid of 1st row, i.e headers.
  getline(reader, line);
  string title, peptide;
  int i = 0;
  while (getline(reader, line)) {
    int pos = line.find('\t');
    title = line.substr(0, pos);
    peptide = line.substr(pos+1);
    
    (*map_title_to_peptide)[title] = peptide;
    (*map_peptide_to_titles)[peptide].emplace_back(title);
  }
  reader.close();
}

void Commons::SpectraSearchBruteForce(
    const HyperParams& params,
    const vector<Spectrum*>& lib_spectra,
    const vector<Spectrum*>& exp_spectra,
    unordered_map<int, vector<int>>& map_spectra_by_charge,
    unordered_map<int, vector<int>>& map_spectra_by_mass,
    vector<int>* top_matches,
    vector<float>* top_scores,
    vector<string>* top_peptides,
    vector<string>* top_titles,
    bool match_isotopic_peak) {

  cout << "SpectraSearchBruteForce with min sim: " << params.min_similarity << endl;
  cout << "SpectraSearchBruteForce with fragment peak precision: " << params.precision << endl;

  omp_set_num_threads(params.threads_to_use);
  cout << "Use #threads: " << params.threads_to_use << endl;

  const int size = exp_spectra.size();
  top_matches->clear(), top_matches->resize(size, -1);
  top_scores->clear(), top_scores->resize(size, -1);
  top_peptides->clear(), top_peptides->resize(size, "NA");
  top_titles->clear(), top_titles->resize(size, "NA");


  float omp_time = 0;
  Utility::Timer timer;
  float num_candidtes = 0;
  float num_candidates_passed_filter = 0;
#pragma omp parallel for reduction(+:omp_time) reduction(+:num_candidtes)
  for (unsigned int i = 0; i < exp_spectra.size(); ++i) {
    //int tid = omp_get_thread_num();
    //cout << "index i: " << i << ", from omp thread: " << tid << endl;;
    Utility::Timer inner_timer;
    const Spectrum& target = *exp_spectra[i];
    vector<float> precursor_list = {target._precursor_mz};
    const vector<float> deviations {-params.precursor_mass_tolerance, 0, 
      params.precursor_mass_tolerance};
    if (match_isotopic_peak) {
      precursor_list.emplace_back(target._precursor_mz - 1. / target._charge);
      precursor_list.emplace_back(target._precursor_mz + 1. / target._charge);
    }
    for (auto precursor : precursor_list) {
      int prev_precursor_id = -1;
      for (auto deviation : deviations) {
        int precursor_id = Utility::Commons::MassToIndex(precursor + deviation);
        if (precursor_id == prev_precursor_id) continue;
        prev_precursor_id = precursor_id;
        auto& indices = map_spectra_by_mass[precursor_id];

        float sim = 0;
        num_candidtes += indices.size();
        for (unsigned int j = 0; j < indices.size(); ++j) {
          const Spectrum& candidate = *lib_spectra[indices[j]];
          if (!Utility::Commons::considerChargeNPrecMass(
                  target._charge, precursor,
                  candidate._charge, candidate._precursor_mz,
                  params.precursor_mass_tolerance)) continue;
          ++num_candidates_passed_filter;
          sim = 1 - Core::Distance::cosine(target, candidate, params.precision);
          if (sim > params.min_similarity && sim > (*top_scores)[i]) {
            (*top_matches)[i] = indices[j];
            (*top_scores)[i] = sim;
            (*top_peptides)[i] = candidate._peptide_raw;
            (*top_titles)[i] = candidate._title;
          }
        }
      }
    }
    omp_time += inner_timer.stop();
  }
  cout << "SpectraSearchBruteForce omp(sum of each thread time) takes secs: " 
      << omp_time << ", omp(general) takes secs: " << timer.stop() << endl; 
  cout << "NOT consider precursor mass and charge, total #candidates: " << num_candidtes << ", avg #candidates: " << num_candidtes*1./size << endl;
  cout << "Consider precursor mass and charge, total #candidates: " << num_candidates_passed_filter << ", avg #candidates: " << num_candidates_passed_filter*1./size << endl;
}

void Commons::CountSpectraPerBucket(
    const vector<Spectrum*>& lib_spectra,
    unordered_map<int, pair<int, int>>* um_table,
    int hash_func_num,
    int hash_dimension,
    int threads_to_use,
    float sample_ratio=1.) {

  const int lib_size = lib_spectra.size();

  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<> dis(0, 1.0);

  const auto hash_table = Core::LSH::generateHashTable(hash_func_num, hash_dimension);

  vector<int> lib_hash_keys(lib_size);

  // Focus on one iteration's distribution of bucket size.
  // unordered_map<int, pair<int, int>> um_table;  // pair<true, decoy>
  
  #pragma omp parallel for num_threads(threads_to_use)
  for (unsigned int i = 0; i < lib_size; ++i) {
      const auto& spectrum = *lib_spectra[i];
    if (dis(gen) <= sample_ratio) { 
      lib_hash_keys[i] = Core::LSH::random_projection(spectrum, hash_table);
    } else {
      lib_hash_keys[i] = -1;
    }
  }

  for (unsigned int i = 0; i < lib_size; ++i) {
    if (-1 == lib_hash_keys[i]) continue;  // -1 means not counted towards subsampled
    const auto& spectrum = *lib_spectra[i];
    //bool decoy = spectrum._peptide_raw == "NA";
    if (spectrum.isDecoy()) {
      (*um_table)[lib_hash_keys[i]].second += 1;
    } else {
      (*um_table)[lib_hash_keys[i]].first += 1;
    }
  }

}

double Commons::AnalyzeBucketSizeFreqDist(
    const HyperParams& params,
    const vector<Spectrum*>& lib_spectra) {
  Utility::Timer timer;

  // Aggregate multiple iterations' distribution of bucket size.
  unordered_map<int, pair<int, int>> um_table;  // <hash_key, pair<#target, #decoy>>

  for (unsigned int it = 0; it < params.iteration; ++it) {
    cout << "assigning spectra to buckets in iteration: " << it << endl;

    CountSpectraPerBucket(lib_spectra,
                          &um_table,
                          params.hash_func_num,
                          params.hash_dimension,
                          params.threads_to_use);

    // Focus on one iteration's distribution of bucket size.
    /*
    string path = params.output_dir.empty() ? "./" : params.output_dir;
    path += "bucket_size_for_h" + to_string(params.hash_func_num) + "i" + to_string(it) + ".csv";
    ofstream writer(path);
    writer << "key,target,decoy" << endl;
    for (const auto& kv : um_table) {
      writer << kv.first << "," << kv.second.first << "," << kv.second.second << endl; 
    }
    writer.close();
    */
  }
  // Aggregate multiple iterations' distribution of bucket size.
  string path = params.output_dir.empty() ? "./" : params.output_dir;
  path += "bucket_size_for_h" + to_string(params.hash_func_num) + "i" + to_string(params.iteration) + ".csv";
  cout << "write bucket size frequency distribution to: " << path << endl;
  ofstream writer(path);
  writer << "key,target,decoy" << endl;
  for (const auto& kv : um_table) {
    writer << kv.first << "," << kv.second.first << "," << kv.second.second << endl; 
  }
  writer.close();
  return timer.stop();
}

double Commons::AnalyzeBucketSizeFreqDistSubsampledSpace(
    const HyperParams& params,
    const vector<Spectrum*>& lib_spectra,
    float sample_ratio=1.) {

  Utility::Timer timer;

  cout << "sample ratio: " << sample_ratio << endl;

  // Aggregate multiple iterations' distribution of bucket size.
  unordered_map<int, pair<int, int>> um_table;  // <hash_key, pair<#target, #decoy>>

  cout << "assigning spectra to buckets" << endl;

  CountSpectraPerBucket(lib_spectra,
                        &um_table,
                        params.hash_func_num,
                        params.hash_dimension,
                        params.threads_to_use,
                        sample_ratio);

  string path = params.output_dir.empty() ? "./" : params.output_dir;
  path += "bucket_size_for_h" + to_string(params.hash_func_num) + ".i1.subsample_ratio" + to_string(sample_ratio) + ".csv";
  ofstream writer(path);
  writer << "key,target,decoy" << endl;
  for (const auto& kv : um_table) {
    writer << kv.first << "," << kv.second.first << "," << kv.second.second << endl; 
  }
  writer.close();
  return timer.stop();
}

double Commons::ApplySingleLSH(
      const HyperParams& params,
      const vector<Spectrum*>& lib_spectra,
      const vector<Spectrum*>& exp_spectra,
      unordered_map<int, unordered_map<int, lsh_table>>* um_lib_hash_table_on_mass_bin,
      unordered_map<int, vector<int>>* um_exp_hash_keys) {

  Utility::Timer timer;
  const int lib_size = lib_spectra.size();
  const int exp_size = exp_spectra.size();

  for (unsigned int it = 0; it < params.iteration; ++it) {
    const auto hash_table = 
        Core::LSH::generateHashTable(params.hash_func_num, params.hash_dimension);
    vector<int> lib_hash_keys(lib_size);
    
    #pragma omp parallel for num_threads(params.threads_to_use)
    for (unsigned int i = 0; i < lib_size; ++i) {
      const auto& spectrum = *lib_spectra[i];
      lib_hash_keys[i] = Core::LSH::random_projection(spectrum, hash_table);
    }
    
    unordered_map<int, lsh_table> local_lsh_table_on_mass_bin;
    for (unsigned int i = 0; i < lib_size; ++i) {
      const auto& spectrum = *lib_spectra[i];
      int bin = Utility::Commons::MassToIndex(spectrum._precursor_mz);
      local_lsh_table_on_mass_bin[bin][lib_hash_keys[i]].emplace_back(i);
    }
    (*um_lib_hash_table_on_mass_bin)[it] = move(local_lsh_table_on_mass_bin);

    vector<int> local_exp_hash_keys(exp_size);
    #pragma omp parallel for num_threads(params.threads_to_use)
    for (unsigned int i = 0; i < exp_size; ++i) {
      const auto& spectrum = *exp_spectra[i];
      local_exp_hash_keys[i] = Core::LSH::random_projection(spectrum, hash_table);
    }
    (*um_exp_hash_keys)[it] = move(local_exp_hash_keys);
  }

  return timer.stop();
}

double Commons::FilterCandidates(
    const HyperParams& params,
    const vector<Spectrum*>& lib_spectra,
    const vector<int>& lib_charge,
    const vector<float>& lib_precursor_mz,
    const vector<Spectrum*>& exp_spectra,
    const unordered_map<int, unordered_map<int, lsh_table>> &um_lib_hash_table_on_mass_bin,
    const unordered_map<int, vector<int>> &um_exp_hash_keys,
    vector<int>* top_matches,
    vector<float>* top_scores,
    vector<string>* top_peptides,
    vector<string>* top_titles,
    bool match_isotopic_peak) {

  omp_set_num_threads(params.threads_to_use);
  cout << "Use #threads: " << params.threads_to_use << endl;

  const int size = exp_spectra.size();

  Utility::CPUTimer cpu_timer; 
  Utility::Timer timer;

  double omp_time = 0;
  float num_candidtes = 0;
  float num_candidates_passed_filter = 0;

  for (int j = 0; j < params.iteration; ++j) {
    const auto& lib_hash_table_on_mass_bin = um_lib_hash_table_on_mass_bin.at(j);
    const auto& exp_hash_keys = um_exp_hash_keys.at(j);
    Utility::Timer inner_timer;

    float num_candidtes_local = 0;

    #pragma omp parallel for reduction(+:num_candidtes) reduction(+:num_candidtes_local)
    for (unsigned int i = 0; i < size; ++i) {
      const Spectrum& target = *exp_spectra[i];
      const int target_charge = target._charge;
      
      vector<float> precursor_list = {target._precursor_mz};
      vector<float> deviations {-params.precursor_mass_tolerance, 0, 
        params.precursor_mass_tolerance};
      
      if (match_isotopic_peak) {
        precursor_list.emplace_back(target._precursor_mz - 1. / target_charge);
        precursor_list.emplace_back(target._precursor_mz + 1. / target_charge);
      }

      for (auto precursor : precursor_list) {
        int prev_precursor_id = -1;
        for (auto deviation : deviations) {
          int precursor_id = Utility::Commons::MassToIndex(precursor + deviation);
          if (precursor_id == prev_precursor_id) continue;
          prev_precursor_id = precursor_id;
        
          const auto& it_precursor = lib_hash_table_on_mass_bin.find(precursor_id);
          if (it_precursor == lib_hash_table_on_mass_bin.end()) continue;

          const int probe_key = exp_hash_keys[i];
          const auto& it_probe = it_precursor->second.find(probe_key);
          if (it_probe == it_precursor->second.end()) continue;

          const auto& probe_values = it_probe->second;

          float sim = 0;
          num_candidtes += probe_values.size();

          num_candidtes_local += probe_values.size();

          for (const auto index : probe_values) {
            const auto& cur = *lib_spectra[index];
            if (!Utility::Commons::considerChargeNPrecMass(
                    target._charge,
                    precursor,
                    //lib_charge[index],
                    cur._charge,
                    //lib_precursor_mz[index],
                    cur._precursor_mz,
                    params.precursor_mass_tolerance)) continue;
            ++num_candidates_passed_filter;

            //const auto& cur = *lib_spectra[index];
            sim = 1 - Core::Distance::cosine(target, cur, params.precision);
            if (sim > params.min_similarity && sim > (*top_scores)[i]) {
              (*top_matches)[i] = index;
              (*top_scores)[i] = sim;
              (*top_peptides)[i] = cur._peptide_raw;
              (*top_titles)[i] = cur._title;
            }
          }
        }
      }
    }
    cout << "[iteration " << j << "] NOT consider precursor mass and charge, total #candidates: " << num_candidtes_local << ", avg #candidates: " << num_candidtes_local*1./size << endl;
    cout << "Time: " << inner_timer.passed() << endl;
    omp_time += inner_timer.stop();
  }

  cout << "FilterCandidates omp(sum of each thread time) takes secs: " 
      << omp_time << ", omp(general) takes secs: " << timer.stop() << endl; 
  cout << "Cpu time: " << cpu_timer.stop() << endl;

  cout << "NOT consider precursor mass and charge, total #candidates: " << num_candidtes << ", avg #candidates: " << num_candidtes*1./size << endl;
  cout << "Consider precursor mass and charge, total #candidates: " << num_candidates_passed_filter << ", avg #candidates: " << num_candidates_passed_filter*1./size << endl;

  return omp_time;
}

void Commons::SpectraSearchLSH(
    const HyperParams& params,
    const vector<Spectrum*>& lib_spectra,
    const vector<int>& lib_charge,
    const vector<float>& lib_precursor_mz,
    const vector<Spectrum*>& exp_spectra,
    vector<int>* top_matches,
    vector<float>* top_scores,
    vector<string>* top_peptides,
    vector<string>* top_titles,
    bool match_isotopic_peak) {
  cout << "SpectraSearchLSH with min sim: " << params.min_similarity << endl;
  cout << "SpectraSearchLSH with fragment peak precision: " << params.precision << endl;
  cout << "SpectraSearchLSH with iteration: " << params.iteration << endl;

  const int size = exp_spectra.size();
  vector<int> local_top_matches(size, -1);
  vector<float> local_top_scores(size, -1);
  vector<string> local_top_peptides(size, "NA");
  vector<string> local_top_titles(size, "NA");

  float time = 0;
  Utility::Timer timer;
  float time_hash = 0, time_filter = 0;
  unordered_map<int, unordered_map<int, lsh_table>> um_lib_hash_table_on_mass_bin;
  unordered_map<int, vector<int>> um_exp_hash_keys;

  cout << "Hashing all spectra with LSH." << endl;
  time_hash = ApplySingleLSH(
      params, 
      lib_spectra, 
      exp_spectra, 
      &um_lib_hash_table_on_mass_bin,
      &um_exp_hash_keys);
  cout << "Hashed all spectra with secs: " << time_hash << endl;

  cout << "Searching nearest neighbor with LSH." << endl;
  time_filter = FilterCandidates(
        params,
        lib_spectra,
        lib_charge,
        lib_precursor_mz,
        exp_spectra,
        um_lib_hash_table_on_mass_bin,
        um_exp_hash_keys,
        &local_top_matches,
        &local_top_scores,
        &local_top_peptides,
        &local_top_titles,
        match_isotopic_peak);
  cout << "Searched nearest neighbor with secs: " << time_filter << endl;

  swap(local_top_matches, *top_matches);
  swap(local_top_scores, *top_scores);
  swap(local_top_peptides, *top_peptides);
  swap(local_top_titles, *top_titles);

  cout << "========= Summary of SpectraSearchLSH =========" << endl;
  cout << "SpectraSearchLSH ApplySingleLSH takes secs: " << time_hash << endl;
  cout << "SpectraSearchLSH FilterCandidates takes secs: " << time_filter << endl;
  cout << "SpectraSearchLSH takes secs: " << timer.stop() << endl;
}

}  // namespace Search
