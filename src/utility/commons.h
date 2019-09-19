#ifndef __UTILITY_COMMONS_H__
#define __UTILITY_COMMONS_H__

#include <cmath>
#include "params.h"

using namespace std;

namespace Utility {
class Commons {
 public:
  static string ReplaceIWithL(const string& peptide);
  static string ReplacePTMOnMWithText(const string& peptide, string text);
  static string ReplacePTMOnCWithText(const string& peptide, string text);
  static int MassToIndex(float mass, float bin_size=1.);
  static inline bool considerChargeNPrecMass(int charge_first, 
                                             float precursor_first, 
                                             int charge_second, 
                                             float precursor_second,
                                             float epsilon) {
    return charge_first == charge_second &&
        fabs(precursor_first - precursor_second) < epsilon;
  }
};

}  // namespace Utility
#endif
