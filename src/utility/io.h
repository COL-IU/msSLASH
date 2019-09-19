#ifndef __UTILITY_IO_H__
#define __UTILITY_IO_H__

#include <algorithm>
#include <chrono>
#include <cassert>
#include <cmath>
#include <ctime>
#include <fcntl.h>
#include <fstream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <unistd.h>
#include <utility>
#include <vector>

#include "../class/spectrum.h"
#include "../class/peak.h"
#include "commons.h"
#include "params.h"

using namespace std;
using namespace Core;

namespace Utility {
enum class PEAK_INTENSITY_RESCALE_METHOD {NONE=0, LOG, SQRT};
class IO {
 public:
  // Use customized char* to float converter, it can support up to 20 decimal
  // places. Changes to [strtof] might be needed to handle higher precision.
  static void ReadSpectraFromMGF(vector<Spectrum*>* indexed_spectra,
      unordered_map<int, vector<int>>* map_spectra_by_charge,
      unordered_map<int, vector<int>>* map_spectra_by_precursor_mass,
      unordered_map<string, vector<int>>* map_peptide_to_indices,
      unordered_map<string, int>* map_ms_title_to_index,
      int* spectra_size, string file_name, float scale, float min_mz,
      float max_mz, float precision, int select_topk, int window_mz,
      PEAK_INTENSITY_RESCALE_METHOD p_inten_rescale_method,
      bool remove_precursor, bool peptide_i2l, bool peptide_ptm_replace, 
      bool filter_unfragmented_ms2, bool verbose=false);


  static void Reset(double* sum_intensity, float* intensity,
                    float* precursor_mz, float* strongest_intensity,
                    int* charge, int* i_peak,
                    string* file_name, string* peptide, string* protein,
                    string* title, string na,
                    Peaks* filtered_peaks);


  // [Deprecated].
  static void ProcessSpectra(vector<Spectrum>* spectra, int start, int end,
      float scale, float min_mz, float precision, int select_topk,
      float window_mz);


  // [Deprecated].
  static void ProcessSpectra(const vector<vector<string>>& chunks, float scale,
      float min_mz, float precision, int select_topk, float window_mz);


  static void SetSpectrum(Spectrum* spectrum,
      bool is_clustered, bool is_consensus, int charge,
      int count, const EmbededPeaks& embeded_peaks, const Peaks& raw_peaks,
      const Peaks& filtered_peaks, string file_name, string protein,
      string raw_peptide, string modified_peptide,  float pre_mz, string title, 
      string component_titles, const vector<float>& top_peak_mz);


  // Normalize peak intensity with the maximum set to 'scale'.
  static void Normalize(Peaks* peaks, float scale);


  // Remove adjacent peaks within mz tolerance, keep the strongest peak.
  static void RemoveAdjacentPeaks(Peaks* peaks, float mz_tolerance);


  // Embed peaks.
  static void Embed(EmbededPeaks* embeded_peaks, const Peaks& peaks,
      float min_mz, float precision, float scale);

  
  // Set consensus spectrum, using just two spectra.
  static void SetConsensus(Spectrum* consensus, const Spectrum& s1,
      const Spectrum& s2, float precision, int topK, float bin_size,
      float min_mz, float scale, string title, string component_titles);


  static void MergeTwoPeaks(const Peaks& p1, const Peaks& p2, Peaks* peaks);


  // Adjust intensity according to frequency.
  static void AdaptPeakIntensities(Peaks* peaks, int nSpectra);


  static void BinTopKPeak(Peaks* top_peaks, const Peaks& peaks, int peaks_size,
      int topK, float bin_size);


  // [Deprecated].
  //static vector<float> SelectTopPeakMZ(
  //    const EmbededPeaks& peaks, float precision,int topK = 5);

  static vector<float> SelectTopPeakMZ(const Peaks& peaks, int topK = 5);
  

  // Customized char* to float converter, up to 20 decimal places.
  static float convert(char const* source, char ** endPtr);
};

}  // namespace Utility
#endif
