#include "io.h"
namespace Utility {

static double const fracMult[] = { 0.0f, 1e-1f, 1e-2f, 1e-3f, 1e-4f, 1e-5f,
    1e-6f, 1e-7f, 1e-8f, 1e-9f, 1e-10f, 1e-11f, 1e-12f, 1e-13f, 1e-14f, 1e-15f,
    1e-16f, 1e-17f, 1e-18f, 1e-19f, 1e-20f};

inline float IO::convert(char const* source, char ** endPtr ) {
  char* end;
  // Change from int to long long, as need large representation on HeLa dataset.
  long long left = strtoll( source, &end, 10 );
  float results = left;
  if ( *end == '.' ) {
      char* start = end + 1;
      long long right = strtoll( start, &end, 10 );
      results += right * fracMult[ end - start ];
  }
  if ( endPtr != nullptr ) {
      *endPtr = end;
  }
  return results;
}

void IO::ReadSpectraFromMGF(vector<Spectrum*>* indexed_spectra,
    unordered_map<int, vector<int>>* map_spectra_by_charge,
    unordered_map<int, vector<int>>* map_spectra_by_precursor_mass,
    unordered_map<string, vector<int>>* map_peptide_to_indices,
    unordered_map<string, int>* map_ms_title_to_index,
    int* spectra_size, string file, float scale, float min_mz,
    float max_mz, float precision, int topK, int bin_size,
    PEAK_INTENSITY_RESCALE_METHOD p_inten_rescale_method, bool remove_precursor, bool remove_precursor_isotopic_peaks, 
    bool peptide_i2l, bool peptide_ptm_replace,
    bool filter_unfragmented_ms2, bool verbose) {

  assert(PEAK_INTENSITY_RESCALE_METHOD::NONE <= p_inten_rescale_method &&
         PEAK_INTENSITY_RESCALE_METHOD::SQRT >= p_inten_rescale_method);

#ifdef FAST_PARSE
  cout << "[FAST_PARSE] Parsing mz and intensity w/ customized parser." << endl;
  cout << "In-house char* to float parser supports up to 20 decimal places." << endl;
  cout << "If this does not suit your needs, compile w/o -D FAST_PARSE" << endl;
#else
  cout << "[STAND_PARSE] Parsing mz and intensity w/ standard parsers." << endl;
#endif
  cout << "Loading spectra from file: " << file << endl;

  auto start_time = chrono::high_resolution_clock::now();

  ifstream reader(file);
  if (!reader.is_open()) {
    cout << "Error! file not open!" << endl;
    exit(-1);
  }

  // Index shall follow next MS/MS spectrum already stored in indexed_spectra.
  const int spectra_cnt_stored = (*indexed_spectra).size();
  int spectra_cnt = 0;

  // Initialize all params.
  Peaks raw_peaks (MAX_PEAK_SIZE);  // Pre-allocate memory.
  Peaks filtered_peaks;  // Filtered peaks using bins.

  string line, tmp;

  const string NA = "NA";
  double sum_intensity = 0;
  float mz, intensity, precursor_mz, strongest_intensity;
  int i_peak, charge;
  string file_name, peptide, protein, title;

  //TODO(leiwang): 50ppm for precursor peak.
  // [abs(threotical - observed) / threotical] * 10^6.
  const float ppm = 50;
  float precursor_peak_tol = 0;
  float precursor_no_water_peak_tol = 0;
  float peak_mz_remove_without_water = -1;


  Reset(&sum_intensity, &intensity,
        &precursor_mz, &strongest_intensity,
        &charge, &i_peak,
        &file_name, &peptide, &protein,
        &title, NA, &filtered_peaks);

  while (getline(reader, line)) {

    // Handle text.
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (0 == line.find(PEPMASS_MGF)) {
      precursor_mz = stof(line.substr(PEPMASS_MGF.length() + 1));
      continue;
    }
    if (0 == line.find(CHARGE_MGF)) {
      charge = stoi(line.substr(CHARGE_MGF.length() + 1));
      continue;
    }
    if (0 == line.find(FILENAME_MGF)) {
      file_name = line.substr(FILENAME_MGF.length() + 1);
      // Get rid of '^M', carriage-and-return created by Windows platform.
      if (file_name.back() == '\r') {
        file_name.pop_back();
      }
      continue;
    }
    if (0 == line.find(SEQ_MGF)) {
      peptide = line.substr(SEQ_MGF.length() + 1);
      // Get rid of '^M', carriage-and-return created by Windows platform.
      if (peptide.back() == '\r') {
        peptide.pop_back();
      }
      continue;
    }
    if (0 == line.find(PEPTIDE_MGF)) {
      peptide = line.substr(PEPTIDE_MGF.length() + 1);
      // Get rid of '^M', carriage-and-return created by Windows platform.
      if (peptide.back() == '\r') {
        peptide.pop_back();
      }
      continue;
    }
    if (0 == line.find(TITLE_MGF)) {
      title = line.substr(TITLE_MGF.length() + 1);
      // Get rid of '^M', carriage-and-return created by Windows platform.
      if (title.back() == '\r') {
        title.pop_back();
      }
      continue;
    }

    if (line[0] < '0' || line[0] > '9') {
      continue;
    }

    peak_mz_remove_without_water =
      precursor_mz - 18./(charge == 0 ? 1 : charge);

    precursor_peak_tol =  precursor_mz * ppm / 1000000;
    precursor_no_water_peak_tol = peak_mz_remove_without_water * ppm / 1000000;

    std::vector<float> precursor_mzs = {precursor_mz};
    if (remove_precursor_isotopic_peaks) {
      precursor_mzs.push_back(precursor_mz + 1. / charge);
      precursor_mzs.push_back(precursor_mz - 1. / charge);
      precursor_mzs.push_back(precursor_mz + 2. / charge);
      precursor_mzs.push_back(precursor_mz - 2. / charge);
    }

    // Read peaks.
    do {
      if (string::npos != line.find(END_IONS_MGF)) {
        break;
      }
      if (line.empty()) {
        continue;
      }

#ifdef FAST_PARSE
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      char* end = nullptr;
      mz = convert(line.c_str(), &end);
      intensity = convert(end, &end);
#else
      char *end = nullptr;
      mz = strtof(line.c_str(), &end);
      intensity = strtof(end, nullptr);
#endif
      if (filter_unfragmented_ms2) {
        strongest_intensity = max(strongest_intensity, intensity);
        sum_intensity += intensity;
      }
      //cout << "[INFO] " << mz << " " << intensity << endl;

      // Rescale intensity with hope to de-emphasize domiant peaks.
      if (PEAK_INTENSITY_RESCALE_METHOD::NONE == p_inten_rescale_method) {
        ;
      } else if (PEAK_INTENSITY_RESCALE_METHOD::LOG == p_inten_rescale_method) {
        intensity = intensity == 0 ? 0 : log(intensity);
      } else if (PEAK_INTENSITY_RESCALE_METHOD::SQRT == p_inten_rescale_method) {
        intensity = intensity <= 0 ? 0 : sqrt(intensity);
      }

      // Filter out peaks out of mz range.
      if (mz < min_mz || mz > max_mz || intensity <= 0) {
        continue;
      }

      // Filter out peaks around precursor_mz and precursor_mz- H2O/charge.
      if (remove_precursor) {
        bool is_precursor_peak = false;
        if (fabs(mz - peak_mz_remove_without_water) < precursor_no_water_peak_tol) {
          is_precursor_peak = true;
        }
        for (const auto& target : precursor_mzs) {
          if (fabs(mz - target) < precursor_peak_tol) {
            is_precursor_peak = true;
            break;
          }
        }
        if (is_precursor_peak) {
          continue;
        }
      }

      // raw_peaks.push_back(Peak(mz, intensity));
      raw_peaks[i_peak++] = Peak(mz, intensity);

    } while(getline(reader, line));

    //TODO(Lei Wang): Adjust the threshold. Maybe 0.9, 0.95.
    if (filter_unfragmented_ms2 &&
        static_cast<double>(strongest_intensity) > sum_intensity * 0.9) {
      continue;
    }

    // This is to handle peaks are not sorted within each spectrum.
    std::sort(raw_peaks.begin(), raw_peaks.begin() + i_peak, [](const Core::Peak& a, const Core::Peak& b) {return a._mz < b._mz;});

    BinTopKPeak(&filtered_peaks, raw_peaks, i_peak, topK, bin_size);
    RemoveAdjacentPeaks(&filtered_peaks, precision);
    Normalize(&filtered_peaks, scale);
    EmbededPeaks embeded_peaks;
    Embed(&embeded_peaks, filtered_peaks, min_mz, precision, scale);

    // Disabled the top_peak features on May 18, 2019.
    //vector<float> top_peak_mz = SelectTopPeakMZ(filtered_peaks);
    vector<float> top_peak_mz;

    // For memory efficient. Discard these 2 lines if memory is sufficient.
    // If descarded, have to clear these vectors before @51L-52L.
    // raw_peaks.clear();
    // filtered_peaks.clear();

    Spectrum* spectrum = new Spectrum();

    string modified_peptide = peptide;
    if (peptide_i2l) modified_peptide = Commons::ReplaceIWithL(modified_peptide);

    // Changes on May 12, 2019 by Lei Wang.
    // If enabled, C+57.021 -> C(C), M+15.995 -> M(O)
    if (peptide_ptm_replace) {
      modified_peptide = Commons::ReplacePTMOnCWithText(modified_peptide, "(C)");
      modified_peptide = Commons::ReplacePTMOnMWithText(modified_peptide, "(O)");
    }

    SetSpectrum(spectrum, false, false, charge, 1, embeded_peaks, Peaks(),
        filtered_peaks, file_name, protein, peptide, modified_peptide,
        precursor_mz, title, title, top_peak_mz);

    (*indexed_spectra).emplace_back(spectrum);

    const int cur_cnt = spectra_cnt_stored + spectra_cnt;

    (*map_spectra_by_charge)[charge].emplace_back(cur_cnt);

    (*map_spectra_by_precursor_mass)[Commons::MassToIndex(precursor_mz)].
        emplace_back(cur_cnt);
    (*map_ms_title_to_index)[title] = cur_cnt;

    if (!peptide.empty()) {
      (*map_peptide_to_indices)[peptide].emplace_back(cur_cnt);
    }

    ++spectra_cnt;
    if (verbose && 0 == spectra_cnt % 100000) {
      cout << "read #spectra: " << spectra_cnt << endl;
    }

    Reset(&sum_intensity, &intensity,
          &precursor_mz, &strongest_intensity,
          &charge, &i_peak,
          &file_name, &peptide, &protein,
          &title, NA,
          &filtered_peaks);
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>
      (end_time - start_time).count();
  cout << "Finish reading " << spectra_cnt << " spectra from file takes secs: " << elapsed_read << endl;

  (*spectra_size) = spectra_cnt;
  reader.close();
}

void IO::Reset(double* sum_intensity, float* intensity,
               float* precursor_mz, float* strongest_intensity,
               int* charge, int* i_peak,
               string* file_name, string* peptide, string* protein,
               string* title, string na, Peaks* filtered_peaks) {
  *sum_intensity = 0.;

  *intensity = 0.;
  *precursor_mz = 0.;
  *strongest_intensity = 0.;

  *charge = 0;
  *i_peak = 0;

  *file_name = na;
  *peptide = na;
  *protein = na;
  *title = na;
  (*filtered_peaks).clear();
}

vector<float> IO::SelectTopPeakMZ(const Peaks& _peaks, int topK) {
  Peaks peaks(_peaks);
  auto cmp = [](const Peak& p1,const Peak& p2)
  {return p1._intensity > p2._intensity;};

  topK = min(topK, (int)peaks.size());
  nth_element(peaks.begin(), peaks.begin() + topK - 1, peaks.end(), cmp);
  vector<float> top_peak_mz(topK);
  for (int i = 0; i < topK; ++i) {
    top_peak_mz[i] = peaks[i]._mz;
  }
  // Sort according to mz ascendingly.
  sort(top_peak_mz.begin(), top_peak_mz.end());
  return top_peak_mz;
}

void IO::SetSpectrum(Spectrum* spectrum,
    bool is_clustered, bool is_consensus, int charge, int count,
    const EmbededPeaks& embeded_peaks, const Peaks& raw_peaks,
    const Peaks& filtered_peaks, string file_name, string protein,
    string raw_peptide, string modified_peptide, float pre_mz,
    string title, string component_titles, const vector<float>& top_peak_mz) {
  (*spectrum)._charge = charge;
  (*spectrum)._component_titles = component_titles;
  (*spectrum)._count = count;
  (*spectrum)._embeded_peaks = embeded_peaks;
  (*spectrum)._filtered_peaks = filtered_peaks;
  (*spectrum)._file_name = file_name;
  (*spectrum)._is_clustered = is_clustered;
  (*spectrum)._is_consensus = is_consensus;
  (*spectrum)._peptide_raw = raw_peptide;
  (*spectrum)._peptide_modified = modified_peptide;
  (*spectrum)._protein = protein;
  (*spectrum)._precursor_mz = pre_mz;
  (*spectrum)._raw_peaks = raw_peaks;
  (*spectrum)._title = title;
  (*spectrum)._top_peak_mz = top_peak_mz;
}

void IO::RemoveAdjacentPeaks(Peaks* peaks, float mz_tolerance) {
  if ((*peaks).empty()) {
    return;
  }
  int i_read, i_write;
  for (i_write = 0, i_read = 1; i_read < (int)(*peaks).size(); ++i_read) {
    auto r_spectrum = (*peaks)[i_read];
    auto& w_spectrum = (*peaks)[i_write];
    if (fabs(r_spectrum._mz - w_spectrum._mz) <= mz_tolerance) {
      if (r_spectrum._intensity > w_spectrum._intensity) {
        w_spectrum = r_spectrum;
      }
    } else {
      (*peaks)[++i_write] = r_spectrum;
    }
  }
  (*peaks).resize(i_write + 1);
}

void IO::Normalize(Peaks* peaks, float scale) {
  float max_intensity = 0;
  for (const auto& peak : *peaks) {
    max_intensity = max(max_intensity, peak._intensity);
  }
  for (auto& peak : *peaks) {
    peak._intensity = peak._intensity / max_intensity * scale;
  }
}

void IO::Embed(EmbededPeaks* embeded_peaks, const Peaks& peaks,
    float min_mz, float precision, float scale) {
  (*embeded_peaks).clear();
  for (const auto& peak : peaks) {
    float mz = peak._mz;
    float intensity = peak._intensity;
    int idx = round((mz - min_mz) / precision);
    (*embeded_peaks).push_back(EmbededPeak(idx, intensity, peak._count));
  }
}

void IO::MergeTwoPeaks (const Peaks& peaks1, const Peaks& peaks2, Peaks* peaks) {

  int i = 0, j = 0;
  int i_len = peaks1.size(), j_len = peaks2.size();
  Peaks local_peaks;
  local_peaks.reserve(i_len + j_len);

  while (i < i_len && j < j_len) {
    if (peaks1[i]._mz == peaks2[j]._mz) {
      //TODO: Need once or twice?
      local_peaks.push_back(peaks1[i]);
      local_peaks.push_back(peaks2[j]);
      ++i;
      ++j;
    } else if (peaks1[i]._mz < peaks2[j]._mz) {
      local_peaks.push_back(peaks1[i]);
      ++i;
    } else {
      local_peaks.push_back(peaks2[j]);
      ++j;
    }
  }

  while (i < i_len) {
    local_peaks.push_back(peaks1[i]);
    ++i;
  }

  while (j < j_len) {
    local_peaks.push_back(peaks2[j]);
    ++j;
  }

  //*peaks = move(local_peaks);
  swap(*peaks, local_peaks);
}

// Merge spectra to build a consensus spectrum.
void IO::SetConsensus(Spectrum* consensus,
    const Spectrum& s1, const Spectrum& s2,
    float precision, int topK, float bin_size, float min_mz, float scale,
    string title, string component_titles) {

  Peaks all_peaks, peaks1, peaks2;
  peaks1.reserve(s1._filtered_peaks.size());
  peaks2.reserve(s1._filtered_peaks.size());

  Peak tmp;
  int count1 = s1._count;
  for (const auto& peak : s1._filtered_peaks) {
    tmp = peak;
    tmp._intensity *= count1;
    peaks1.push_back(tmp);
  }

  int count2= s2._count;
  for (const auto& peak : s2._filtered_peaks) {
    tmp = peak;
    tmp._intensity *= count2;
    peaks2.push_back(tmp);
  }

  MergeTwoPeaks(peaks1, peaks2, &all_peaks);

  int total_count = count1 + count2;
  float ave_precursor_mz =
    (s1._precursor_mz * count1 + s2._precursor_mz * count2 ) / total_count;

  int size = s1._filtered_peaks.size() + s2._filtered_peaks.size();

  //TODO: Better options? Can change to while loop for more iterations.
  Peaks new_peaks(size);
  float threshold = precision;
  new_peaks[0] = all_peaks[0];
  int i_last = 0;
  for (int i = 1; i < all_peaks.size(); ++i) {
    if (all_peaks[i]._mz - new_peaks[i_last]._mz <= threshold) {
      auto& last_peak = new_peaks[i_last];
      const auto& to_merge_peak = all_peaks[i];
      float sum_intensity = last_peak._intensity + to_merge_peak._intensity;

      // Set summed inten, weighted mz, count, i.e(spectrum that has it peak)
      last_peak._mz = last_peak._mz * last_peak._intensity / sum_intensity +
        to_merge_peak._mz * to_merge_peak._intensity / sum_intensity;
      last_peak._intensity = sum_intensity;

    } else {
      // new_peaks.push_back(all_peaks[i]);
      new_peaks[++i_last] = all_peaks[i];
    }
  }
  new_peaks.resize(++i_last);

  Peaks filtered_peaks;
  BinTopKPeak(&filtered_peaks, new_peaks, i_last, topK, bin_size);
  RemoveAdjacentPeaks(&filtered_peaks, precision);
  Normalize(&filtered_peaks, scale);
  EmbededPeaks embeded_peaks;
  Embed(&embeded_peaks, filtered_peaks, min_mz, precision, scale);

  vector<float> top_peak_mz = SelectTopPeakMZ(filtered_peaks);

  // For memory efficient. Discard these 2 lines if memory is sufficient.
  new_peaks.clear();
  // filtered_peaks.clear();

  int charge = (0 == s1._charge ? s2._charge: s1._charge);

  string file_name = "NA", protein = "NA", raw_peptide = "NA",
         modified_peptide = "NA";

  SetSpectrum(consensus, false, true, charge, total_count, embeded_peaks,
      new_peaks, filtered_peaks, file_name, protein, raw_peptide,
      modified_peptide, ave_precursor_mz, title, component_titles, top_peak_mz);

}

// Adapt the peak intensities in consensus peak using the following formula:
// I = I * (0.95 + 0.05 * (1 + pi)^5)
void IO::AdaptPeakIntensities(Peaks* peaks, int nSpectra) {
  for (auto& peak : *peaks) {
    float pi = (float)peak._count / nSpectra;
    //TODO: Plan to deprecate it.
    // assert(pi <= 1);
    if (pi > 1) {
      cout << peak << endl;
      cout << "nSpectra: " << nSpectra << endl;
      cout << "Pi greater than 1!" << endl;
      pi = 1.;
    }
    float new_intensity = peak._intensity * (0.95 + 0.05 * pow(1 + pi, 5));
    peak._intensity = new_intensity;
  }
}

void IO::BinTopKPeak(Peaks* top_peaks, const Peaks& peaks, int peaks_size,
    int topK, float bin_size) {
  // Assume that peaks were already sorted ascendingly w.r.t mz.
  float max_mz = peaks[peaks_size - 1]._mz;
  // Sort peaks descendingly w.r.t peak's intensity.
  auto cmp =
    [](const Peak& p1, const Peak& p2){return p1._intensity > p2._intensity;};

  *top_peaks = Peaks((round(max_mz / bin_size) + 1 ) * topK);
  Peaks current_peaks(peaks_size);
  int i_write = 0;

  int last_index = 0;
  float current_max_mz = 0;
  for (float min_mz = 0; min_mz < max_mz; min_mz += bin_size) {
    current_max_mz = min_mz + bin_size;
    int i_current_peak = 0;
    for (; last_index < peaks_size; ++last_index) {
      const auto& peak = peaks.at(last_index);
      if (peak._mz > current_max_mz) {
        break;
      }
      current_peaks[i_current_peak++] = peak;
    }

    int ele_num = min(i_current_peak, topK);
    //TODO: Should be .begin() + ele_num - 1.
    nth_element(current_peaks.begin(), current_peaks.begin() + ele_num - 1,
        current_peaks.begin() + i_current_peak, cmp);
    for (int i = 0; i < ele_num; ++i) {
      (*top_peaks)[i_write++] = current_peaks[i];
    }
  }
  (*top_peaks).resize(i_write);
  sort((*top_peaks).begin(), (*top_peaks).end());
}

}  // namespace Utility
