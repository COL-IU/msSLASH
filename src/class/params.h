#ifndef __PARAMS_H__
#define __PARAMS_H__

#include <vector>
#include <string>

#include "peak.h"
using std::vector;
using std::string;

namespace Core {
typedef vector<EmbededPeak> EmbededPeaks;
typedef vector<Peak> Peaks;

// Constants for spectrum constructor
const string NA = "NA";
const string DECOY = "_DECOY";
}  // namespace Core
#endif
