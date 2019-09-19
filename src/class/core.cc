#include "core.h"
namespace Core {
hashFunction LSH::generateNormalHashFunc(int size) {
  hashFunction function(size, 0);

  random_device rd;
  mt19937 gen(rd());

  normal_distribution<> dis(0, 1);
  for (int i = 0; i < size; ++i) {
    function[i] = dis(gen);
  }
  
  return function;
}

hashFunction LSH::generateUniformHashFunc(int size, float min, float max) {
  hashFunction function(size, 0);

  random_device rd;
  mt19937 gen(rd());

  uniform_real_distribution<> dis(min, max);
  for (int i = 0; i < size; ++i) {
    function[i] = dis(gen);
  }
  
  return function;
}

hashTable LSH::generateHashTable(int num_hash_func, int dim_hash_func) {
  hashTable table;
  for (int j = 0; j < num_hash_func; ++j) {
    table.emplace_back(generateNormalHashFunc(dim_hash_func));
  }
  return table;
}

/*
int LSH::random_projection(const EmbededPeaks& embeded_peaks,
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

int LSH::random_projection(const EmbededPeaks& embeded_peaks,
    const hashTable& table, bool binarized) {
  int key = 0;
  for (auto& function : table) {
    key = key * 2 + random_projection(embeded_peaks, function, binarized);
  }
  return key;
}

int LSH::random_projection(const Spectrum& spectrum,
    const hashTable& table, bool binarized) {
  return random_projection(spectrum._embeded_peaks, table, binarized);
}
*/

//TODO: Test p_stable functions.
string LSH::p_stable(const EmbededPeaks& embeded_peaks,
    const hashFunction& function, float b, float r) {
  float sum = 0;
  for (const auto& peak : embeded_peaks) {
    sum += peak._intensity * function[peak._idx];
  }
  return to_string(int((sum + b) / r));
}
string LSH::p_stable(const EmbededPeaks& embeded_peaks,
    const hashTable& table, float b, float r) {
  string key;
  for (const auto& function : table) {
    key += p_stable(embeded_peaks, function, b, r);
  }
  return key;
}
}

