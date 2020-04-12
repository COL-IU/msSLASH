#ifndef  __SEARCH_H__
#define  __SEARCH_H__
#include <ctime>
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
#include "../class/core.h"
#include "../class/distance.h"
#include "../utility/io.h"
#include "../utility/timer.h"
#include "../utility/commons.h"

using namespace std;

namespace Search {
class Commons {
 public:
  static string getTimeAndDate();

  // I -> L, C+57.021 -> C(C), M+15.995 -> M(O)
  static string AdjustPeptide(const string& peptide);

  static void ReadSpectraHelper(
    const HyperParams& params,
    const string& file, 
    vector<Spectrum*>* spectra, 
    unordered_map<int, vector<int>>* map_spectra_by_charge,
    unordered_map<int, vector<int>>* map_spectra_by_mass,
    unordered_map<string, vector<int>>* map_peptide_to_indices,
    unordered_map<string, int>* map_ms_titles_to_index,
    bool peptide_i2l, bool peptide_ptm_replace);

  static void 
      ReleaseMemoryHelper(string spectra_name, vector<Spectrum*>* spectra);

  static void ReadMsgfIds(
    string fileName, 
    unordered_map<string, string>* map_title_to_peptide, 
    unordered_map<string, vector<string>>* map_peptide_to_titles);

  static void SpectraSearchBruteForce(
    const HyperParams& params,
    const vector<Spectrum*>& lib_spectra,
    const vector<Spectrum*>& exp_spectra,
    unordered_map<int, vector<int>>& map_spectra_by_charge,
    unordered_map<int, vector<int>>& map_spectra_by_mass,
    vector<int>* top_matches,
    vector<float>* top_scores,
    vector<string>* top_peptides,
    vector<string>* top_titles,
    bool match_isotopic_peak);

  static double ApplySingleLSH(
      const HyperParams& params,
      const vector<Spectrum*>& lib_spectra,
      const vector<Spectrum*>& exp_spectra,
      unordered_map<int, unordered_map<int, lsh_table>>* um_lib_hash_table_on_mass_bin,
      unordered_map<int, vector<int>>* um_exp_hash_keys);

   static void CountSpectraPerBucket(
      const vector<Spectrum*>& lib_spectra,
      unordered_map<int, pair<int, int>>* um_table,
      int hash_func_num,
      int hash_dimension,
      int threads_to_use,
      float sample_ratio);

  // Anlyze the frequency distribution of number of spectra inside per bucket
  static double AnalyzeBucketSizeFreqDist(
      const HyperParams& params,
      const vector<Spectrum*>& lib_spectra);

  // Anlyze the frequency distribution of number of spectra inside per bucket
  static double AnalyzeBucketSizeFreqDistSubsampledSpace(
      const HyperParams& params,
      const vector<Spectrum*>& lib_spectra,
      float sample_ratio);

  static double FilterCandidates(
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
    bool match_isotopic_peak);

  static void SpectraSearchLSH(
    const HyperParams& params,
    const vector<Spectrum*>& lib_spectra,
    const vector<int>& lib_charge,
    const vector<float>& lib_precursor_mz,
    const vector<Spectrum*>& exp_spectra,
    vector<int>* top_matches,
    vector<float>* top_scores,
    vector<string>* top_peptides,
    vector<string>* top_titles,
    bool match_isotopic_peak);
};

}  // namespace Search
#endif
