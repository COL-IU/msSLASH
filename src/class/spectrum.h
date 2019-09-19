#ifndef __CORE_SPECTRUM_H__
#define __CORE_SPECTRUM_H__

#include <cmath>
#include <iostream>
#include <string>

#include "params.h"

using std::endl;
using std::fabs;
using std::ostream;
using std::string;

namespace Core {
class Spectrum {
 public:
  bool _is_clustered;
  bool _is_consensus;
  float _precursor_mz;
  int _charge;
  // Number of spectra contributing to this spectrum.
  int _count;
  // Ascendingly sorted mz of topK strongest peaks.
  vector<float> _top_peak_mz;

  // Store the titles of its component spectra.
  // For instance:
  //  1. If this spectrum is not a consensus spectrum,
  //    then _title = _component_titles.
  //  2. If this spectrum is a consensus spectrum composed of s_a, s_b, s_c,
  //    then _component_titles =
  //    ";".join(_title of s_a, _title of s_b, _title of s_c).
  string _component_titles;
  string _file_name;
  string _peptide_raw;  // Extracted directly from mgf file, if any.
  string _peptide_modified;  // Possible modifications might contain I2L, PTMs.
  string _protein;
  string _title;
  // Embedded peaks, with mz-to-idx and intensity.
  EmbededPeaks _embeded_peaks;
  // Filtered peaks, such as Top5 per 100 Da.
  Peaks _filtered_peaks;
  // Raw peaks, with mz in [min_mz, max_mz].
  Peaks _raw_peaks;
  
  Spectrum(): _is_clustered(false), _is_consensus(false), _precursor_mz(0), 
    _charge(0), _count(1), _file_name("NA"), _peptide_raw("NA"), 
    _peptide_modified("NA"), _protein("NA"), _title("NA"), 
    _component_titles("NA") {};

  ~Spectrum() {}

  Spectrum(const Spectrum& sp) {
    _charge = sp._charge;
    _component_titles = sp._component_titles;
    _count = sp._count;
    _embeded_peaks = sp._embeded_peaks;
    _filtered_peaks = sp._filtered_peaks;
    _is_clustered = sp._is_clustered;
    _is_consensus = sp._is_consensus;
    _file_name = sp._file_name;
    _peptide_raw = sp._peptide_raw;
    _peptide_modified = sp._peptide_modified;
    _protein = sp._protein;
    _precursor_mz = sp._precursor_mz;
    _raw_peaks = sp._raw_peaks;
    _title = sp._title;
    _top_peak_mz = sp._top_peak_mz;
  }

  Spectrum& operator=(const Spectrum& sp) {
    _charge = sp._charge;
    _component_titles = sp._component_titles;
    _count = sp._count;
    _embeded_peaks = sp._embeded_peaks;
    _filtered_peaks = sp._filtered_peaks;
    _is_clustered = sp._is_clustered;
    _is_consensus = sp._is_consensus;
    _file_name = sp._file_name;
    _peptide_raw = sp._peptide_raw;
    _peptide_modified = sp._peptide_modified;
    _protein = sp._protein;
    _precursor_mz = sp._precursor_mz;
    _raw_peaks = sp._raw_peaks;
    _title = sp._title;
    _top_peak_mz = sp._top_peak_mz;
    return *this;
  }
  
  friend ostream& operator<<(ostream& os, const Spectrum& spectrum);

  bool shareTopPeaks(const Spectrum& other, float epsilon) {
    int i = 0, j = 0;
    while (i < _top_peak_mz.size() && j < other._top_peak_mz.size()) {
      if (fabs(_top_peak_mz[i] - other._top_peak_mz[j]) < epsilon) {
        return true;
      }
      if (_top_peak_mz[i] < other._top_peak_mz[j]) {
        ++i;
      } else {
        ++j;
      }
    }
    return false;
  }

  // Share similar precursor mass or not.
  bool shareSimPrecMass(const Spectrum& other, float epsilon) const {
    return fabs(_precursor_mz - other._precursor_mz) < epsilon;
  }

  // Share same charge && similar precursor mass or not.
  bool considerChargeNPrecMass(const Spectrum& other, float epsilon) const {
    if (_charge != other._charge) return false;
    return shareSimPrecMass(other, epsilon);
  }

  bool considerChargeNPrecMass(int charge, float precursor_mz, float epsilon) const {
    if (_charge != charge) return false;
    return fabs(_precursor_mz - precursor_mz) < epsilon;
  }
};

inline ostream& operator<<(ostream& os, const Spectrum& spectrum) {
  os << "is_clustered: " << spectrum._is_clustered << endl;
  os << "is_consensus: " << spectrum._is_consensus << endl;
  os << "precursor_mz: " << spectrum._precursor_mz << endl;
  os << "charge: " << spectrum._charge << endl;
  os << "file name: " << spectrum._file_name << endl;
  os << "raw peptide: " << spectrum._peptide_raw << endl;
  os << "modified peptide: " << spectrum._peptide_modified << endl;
  os << "protein: " << spectrum._protein << endl;
  os << "title: " << spectrum._title << endl;
  os << "component titles: " << spectrum._component_titles << endl;
  os << "raw peaks: " << endl;
  for (const auto& peak : spectrum._raw_peaks) {
    os << peak << endl;
  }
  os << endl;
  os << "filtered peaks: " << endl;
  for (const auto& peak : spectrum._filtered_peaks) {
    os << peak << endl;
  }
  os << endl;
  os << "top peaks' mz:" << endl;
  for (auto mz : spectrum._top_peak_mz) {
    os << mz << endl;
  }
  os << endl;
  os << "embeded_peaks: " << endl;
  for (const auto& peak : spectrum._embeded_peaks) {
    os << peak << endl;
  }
  os << endl;
  return os;
}

}  // namespace Core
#endif
