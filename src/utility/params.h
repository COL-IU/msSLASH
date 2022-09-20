#ifndef __UTILITY_PARAMS_H__
#define __UTILITY_PARAMS_H__

#include <utility>
#include <regex>
#include <string>
#include <vector>

using std::pair;
using std::string;
using std::vector;
using std::regex;

namespace Utility {
// Constants for spectrum in .sptxt file.
const string BINARYFILEOFFSET_SPTXT = "BinaryFileOffset";
const string NAME_SPTXT = "Name";
const string NUMPEAKS_SPTXT = "NumPeaks";
const string PRECURSORMZ_SPTXT = "PrecursorMZ";

// Constants for spectrum in .msp file.
const string NAME_MSP = "Name";
const string NUMPEAKS_MSP = "Num peaks";
const string PRECURSORMZ_MSP = "Parent";
// const string BINARYFILEOFFSET_MSP = "BinaryFileOffset";

// Constants for SpectraST search result in .pep.xml file.
const string LIB_FILE_OFFSET_XML = "lib_file_offset";
const string SPECTRUM_XML = "spectrum=\"";
const string VALUE_XML = "value";

// Constants for spectrum in .mgf file.
const string BEGIN_IONS_MGF = "BEGIN IONS";
const string CHARGE_MGF = "CHARGE";
const string END_IONS_MGF = "END IONS";
const string PEPMASS_MGF = "PEPMASS";
const string PEPTIDE_MGF = "PEPTIDE";
const string TITLE_MGF = "TITLE"; 
const string FILENAME_MGF = "FILENAME"; 
const string SEQ_MGF = "SEQ";

// [Deprecated] Constants for modificatons for peptide.
const string OXIDATION_M = "(O)"; 
const string Carbamidomethyl_C = "+57.021";

// Regular expression for the PTMs on M and C, respectively.
// const regex regPTM_M("M(\\(O\\)|\\+?[[:digit:]]+\\.?[[:digit:]]+)");  // M(O) or M+15.995
const regex regPTM_M("M\\+15\\.995");  // M+15.995 for MSGF searching result.
//const regex regPTM_C("C(\\+?[[:digit:]]+\\.?[[:digit:]]+)");  // C+57.021
const regex regPTM_C("C\\+57\\.021");  // C+57.021 for MSGF searching result.

// Constants for raw peaks size for pre-allocating memory.
const int MAX_PEAK_SIZE = 1e5;
}  // namespace Utility
#endif
