#include "commons.h"
namespace Utility {

string Commons::ReplaceIWithL(const string& peptide) {
  string pep = peptide;
  for (size_t i = 0; i < pep.length(); ++i) {
    if ('I' == pep[i]) pep[i] = 'L';
  }
  return pep;
}

string Commons::ReplacePTMOnMWithText(const string& peptide, string text) {
  return regex_replace(peptide, Utility::regPTM_M, "M" + text);
}

string Commons::ReplacePTMOnCWithText(const string& peptide, string text) {
  return regex_replace(peptide, Utility::regPTM_C, "C" + text);
}

int Commons::MassToIndex(float mass, float bin_size) {
  return int(round(mass / bin_size));
}

}  // namespace Utility
