#include <ctime>
#include <cmath>
#include <iostream>
#include <map>
#include <mutex>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <string>
#include <thread>
#include <omp.h>

#include "../class/hyperparams.h"
#include "../class/params.h"
#include "../class/spectrum.h"
#include "../class/core.h"
#include "../class/distance.h"
#include "../utility/io.h"
#include "../utility/commons.h"
#include "../utility/cmdparser.h"
#include "./commons.h"

using namespace std;

void SetParameters(
    const cli::Parser& parser, 
    Core::HyperParams* params, 
    string* lib_file, 
    string* exp_file) {
  (*params).filter_unfragmented_ms2 = parser.get<bool>("u");
  (*params).hash_func_num = parser.get<int>("n");
  (*params).iteration = parser.get<int>("i");
  (*params).peak_intensity_rescale_method = parser.get<int>("r");
  (*params).min_mz = parser.get<float>("a");
  (*params).min_similarity = parser.get<float>("s");
  (*params).threads_to_use = parser.get<int>("t");
  (*params).precision = parser.get<float>("c");
  (*params).precursor_mass_tolerance = parser.get<float>("m");
  (*params).resizeHashDim();  // Fixed on 09062019

  (*params).msSLASH_tsv_file = parser.get<string>("o");
  (*params).use_precision_one_thompson = parser.get<float>("precision_1th");

  *lib_file = parser.get<string>("l");
  *exp_file = parser.get<string>("e");

  if ((*params).msSLASH_tsv_file.empty()) {
    (*params).msSLASH_tsv_file = 
        (*exp_file).substr(0, (*exp_file).rfind(".")) + ".msSLASH.tsv";
  }
}

void configure_parser(cli::Parser& parser) {
  parser.set_optional<bool>("u", "unfragment", false, "[Bool] Filter unfragmented ms2.");
  parser.set_optional<bool>("precision_1th", "precision_1th", false, "use 1th=1/charge for precision");

  parser.set_optional<int>("i", "iteration", 100, "[Int] iteration for searching with LSH.");
  parser.set_optional<int>("n", "hash_func_num", 8, "[Int] hash functions for LSH.");
  parser.set_optional<int>("r", "rescale", 1, "[Int] PEAK_INTENSITY_RESCALE_METHOD.");
  parser.set_optional<int>("t", "threads", 20, "[Int] num of threads to use.");

  parser.set_optional<float>("a", "min_mz", 200, "[Float] min mz for peaks.");
  parser.set_optional<float>("c", "precision", 0.5, "[Float] fragment precision.");
  parser.set_optional<float>("m", "precursor_mass_tolerance", 0.05, "[Float] precursor mass tolerance in Da.");
  parser.set_optional<float>("s", "similarity", 0., "[Float] similarity threshold .");

  parser.set_optional<string>("o", "out_file", "", "[String] out file contanining msSLASH searching results.");

  parser.set_required<string>("d", "decoy", "[String] decoy library mgf file.");
  parser.set_required<string>("e", "experimental", "[String] experimental mgf file.");
  parser.set_required<string>("l", "library", "[String] library mgf file.");
}

int main (int argc, char *argv[]) {
  cli::Parser parser(argc, argv);
  configure_parser(parser);
  parser.run_and_exit_if_error();

  Core::HyperParams params;
  string lib_file, exp_file;
  string decoy_file = parser.get<string>("d");

  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  SetParameters(parser, &params, &lib_file, &exp_file);
  cout << params << endl;
  cout << "target library mgf file: " << lib_file << endl;
  cout << "decoy library mgf file: " << decoy_file<< endl;
  cout << "experimental dataset mgf file: " << exp_file << endl;
  cout << endl;

  bool peptide_i2l = true;  // peptide I2L.
  bool peptide_ptm_replace = true;  // C(C) and M(O)

  bool match_isotopic_peak = true;
  // bool match_isotopic_peak = false;

  // Read target library MGF files.
  vector<Spectrum*> lib_spectra;
  unordered_map<int, vector<int>> map_lib_spectra_by_charge;
  unordered_map<int, vector<int>> map_lib_spectra_by_mass;
  unordered_map<string, vector<int>> map_lib_peptide_to_indices;
  unordered_map<string, int> map_lib_ms_title_to_index;
  Search::Commons::ReadSpectraHelper(params, lib_file, &lib_spectra, 
                                     &map_lib_spectra_by_charge,
                                     &map_lib_spectra_by_mass, 
                                     &map_lib_peptide_to_indices,
                                     &map_lib_ms_title_to_index,
                                     peptide_i2l,
                                     peptide_ptm_replace);

  // Read decoy library MGF files.
  Search::Commons::ReadSpectraHelper(params, decoy_file, &lib_spectra, 
                                     &map_lib_spectra_by_charge,
                                     &map_lib_spectra_by_mass, 
                                     &map_lib_peptide_to_indices,
                                     &map_lib_ms_title_to_index,
                                     peptide_i2l,
                                     peptide_ptm_replace);

  unordered_map<string, vector<int>> map_lib_modified_peptide_to_indices;
  for (const auto& kv : map_lib_peptide_to_indices) {
    string peptide = kv.first;
    if (peptide_i2l && peptide_ptm_replace) {
      peptide = Search::Commons::AdjustPeptide(peptide);
    }
    map_lib_modified_peptide_to_indices[peptide] = kv.second;
  }

  // Read EXPERIMENTAL MGF files.
  vector<Spectrum*> exp_spectra;
  unordered_map<int, vector<int>> map_exp_spectra_by_charge;
  unordered_map<int, vector<int>> map_exp_spectra_by_mass;
  unordered_map<string, vector<int>> map_exp_peptide_to_indices;
  unordered_map<string, int> map_exp_ms_title_to_index;
  Search::Commons::ReadSpectraHelper(params, exp_file, &exp_spectra, 
                                     &map_exp_spectra_by_charge,
                                     &map_exp_spectra_by_mass, 
                                     &map_exp_peptide_to_indices,
                                     &map_exp_ms_title_to_index,
                                     peptide_i2l,
                                     peptide_ptm_replace);

  // Code should be added after this line.
  
  vector<int> top_matches;
  vector<float> top_scores;
  vector<string> top_raw_peptides;
  vector<string> top_titles;
#ifdef __BRUTEFORCE__
  cout << "SpectraSearchBruteForce enabled." << endl; 
  cout << "msSLASH results write to file: "  << params.msSLASH_tsv_file << endl;
  Search::Commons::SpectraSearchBruteForce(
      params, 
      lib_spectra, 
      exp_spectra, 
      map_lib_spectra_by_charge, 
      map_lib_spectra_by_mass, 
      &top_matches,
      &top_scores, 
      &top_raw_peptides, 
      &top_titles,
      match_isotopic_peak);
#else
  cout << "SpectraSearchLSH enabled." << endl; 
  cout << "msSLASH results write to file: "  << params.msSLASH_tsv_file << endl;
  vector<int> lib_charge;
  vector<float> lib_precursor_mz;
  for (const auto& spectrum : lib_spectra) {
    lib_charge.emplace_back(spectrum->_charge);
    lib_precursor_mz.emplace_back(spectrum->_precursor_mz);
  }
  Search::Commons::SpectraSearchLSH(
      params,
      lib_spectra,
      lib_charge,
      lib_precursor_mz,
      exp_spectra,
      &top_matches,
      &top_scores, 
      &top_raw_peptides, 
      &top_titles,
      match_isotopic_peak);
#endif

  ofstream writer(params.msSLASH_tsv_file);
  writer << "Index\t"
      "Title\t"
      "Charge\t"
      "TopMatch\t"
      "TopScore\t"
      "TopPep\t"
      << endl;

  for (unsigned int i = 0; i < top_matches.size(); ++i) {
    const auto& spectrum = *exp_spectra[i];
    writer << i << "\t"
        << spectrum._title << "\t" 
        << spectrum._charge << "\t" 
        << top_matches[i] << "\t" 
        << top_scores[i] << "\t" 
        << top_raw_peptides[i] << "\t"
        << endl;
  }
  writer.close();


  // Code should be added before this line.
  
  // Releasing memory of spectra.
  cout << endl;
  Search::Commons::ReleaseMemoryHelper("lib_spectra", &lib_spectra);
  Search::Commons::ReleaseMemoryHelper("exp_spectra", &exp_spectra);

  return 0;
}

