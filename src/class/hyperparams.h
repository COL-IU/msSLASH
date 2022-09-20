#ifndef __CORE_HYPERPARAMS_H__
#define __CORE_HYPERPARAMS_H__

#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;
using std::ostream;

namespace Core {
class HyperParams{
 public:
  bool remove_precursor;
  bool filter_unfragmented_ms2;
  bool use_precision_one_thompson;

  float precursor_mass_tolerance;
  float max_mz;
  float min_mz;
  float min_similarity;
  float mz_scale;
  float precision;
  float shared_peak_mz_epsilon;
  float union_find_merging_threshold;

  int hash_func_num;
  int hash_dimension;
  int iteration;
  int peak_intensity_rescale_method;
  int select_topk;
  int threads_to_use;
  int window_mz;

  string file_name;
  string result_path;
  string result_prefix;
  string cs_path;
  string msgf_tsv_file;
  string msgf_tsv_psm_file;
  string msSLASH_tsv_file;
  string cluster_file;
  string suffix;
  string output_dir;

  HyperParams() {
    cout << "--Default Constructor for HyperParams called.--" << endl; 
    iteration = 100;
    hash_func_num = 10;
    peak_intensity_rescale_method = 1;
    precursor_mass_tolerance = 0.05;
    max_mz = 2000;
    union_find_merging_threshold = 0.1;
    min_mz = 200;                                                                                                  
    min_similarity = 0.;
    mz_scale = 1000;
    precision = 0.8;  // Fragment precision.
    remove_precursor = true;
    filter_unfragmented_ms2 = false;
    shared_peak_mz_epsilon = 0.2;
    select_topk = 5;
    threads_to_use = 20;
    window_mz = 100;
    hash_dimension = int((max_mz - min_mz) / precision) + 1;
  }
  void resizeHashDim() {
    try{
      hash_dimension = int((max_mz - min_mz) / precision) + 1;
    } catch (const std::exception& e) {
      std::cout << "a standard exception was caught, with message '"
          << e.what() << "'\n";
    }
  }

  friend ostream& operator<<(ostream& os, const HyperParams& params);
};

inline ostream& operator<<(ostream& os, const HyperParams& params) {
  os << "--Output HyperParams Begins--" << endl;

  os << "[bool] remove_precursor: " << params.remove_precursor << endl;
  os << "[bool] filter_unfragmented_ms2: " << params.filter_unfragmented_ms2 << endl;

  os << "[float] precursor_mass_tolerance: " << params.precursor_mass_tolerance << endl;
  os << "[float] max_mz: " << params.max_mz << endl;
  os << "[float] min_mz: " << params.min_mz << endl;
  os << "[float] min_similarity: " << params.min_similarity << endl;
  os << "[float] mz_scale: " << params.mz_scale << endl;
  os << "[float] precision: " << (params.use_precision_one_thompson ? "1Th" : std::to_string(params.precision)) << endl;
  os << "[float] shared_peak_mz_epsilon: " << params.shared_peak_mz_epsilon << endl;
  os << "[float] union_find_merging_threshold: " << params.union_find_merging_threshold << endl;

  os << "[int] hash_func_num: " << params.hash_func_num << endl;
  os << "[int] hash_dimension: " << params.hash_dimension << endl;
  os << "[int] iteration: " << params.iteration << endl;
  os << "[int] peak_intensity_rescale_method (0: NONE; 1:LOG; 2:SQRT): " << params.peak_intensity_rescale_method << endl;
  os << "[int] select_topk: " << params.select_topk << endl;
  os << "[int] threads_to_use: " << params.threads_to_use << endl;
  os << "[int] window_mz: " << params.window_mz << endl;
  os << "[string] output_dir: " << params.output_dir << endl;
  os << "[string] file_name: " << params.file_name << endl;
  os << "[string] result_path: " << params.result_path << endl;
  os << "[string] result_prefix: " << params.result_prefix << endl;
  os << "[string] cs_path: " << params.cs_path << endl;
  os << "[string] msgf_tsv_file: " << params.msgf_tsv_file << endl;
  os << "[string] msgf_tsv_psm_file: " << params.msgf_tsv_psm_file << endl;
  os << "[string] msSLASH_tsv_file: " << params.msSLASH_tsv_file << endl;
  os << "[string] cluster_file: " << params.cluster_file << endl;
  os << "[string] suffix: " << params.suffix << endl;

  os << "--Output HyperParams Ends--" << endl;
  return os;
}

}  // namespace Core
#endif
