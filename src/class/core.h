#ifndef __CORE_H__
#define __CORE_H__

#include <algorithm>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "params.h"
#include "spectrum.h"

using namespace std;

namespace Core {
typedef vector<float> hashFunction;
typedef vector<hashFunction> hashTable;
typedef unordered_map<int, vector<int> > lsh_table;
class LSH {
 public:
  static hashFunction generateNormalHashFunc(int size);
  static hashFunction generateUniformHashFunc(int size, float min, float max);
  static hashTable generateHashTable(int num_hash_func, int dim_hash_func); 
  
  // Different LSHs.
  static inline int random_projection(const EmbededPeaks& embeded_peaks,
      const hashFunction& function, bool binarized) {
    float sum = 0;
    for (auto& peak : embeded_peaks) {
      if (binarized) {
        sum += function[peak._idx];
      } else { 
        sum += peak._intensity * function[peak._idx];
      }
    }
    return sum >=0 ? 1 : 0;
  }

  static inline int random_projection(const EmbededPeaks& embeded_peaks,
      const hashTable& table, bool binarized = false) {
    int key = 0;
    for (auto& function : table) {
      key = key * 2 + random_projection(embeded_peaks, function, binarized);
    }
    return key;
  }

  static inline int random_projection(const Spectrum& spectrum,
      const hashTable& table, bool binarized = false) {
    return random_projection(spectrum._embeded_peaks, table, binarized);
  }

  static string p_stable(const EmbededPeaks& embeded_peaks,
      const hashFunction& function, float b, float r);
  static string p_stable(const EmbededPeaks& embeded_peaks,
      const hashTable& table, float b, float r);
};

}  // namespace Core
#endif
