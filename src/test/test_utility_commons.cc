#include <iostream>
#include "../utility/commons.h"

using namespace std;

int main() {

  auto inputs = {"AVTNHISVYC+57.021STK", "SM(O)TSAHGSASVNS(O)R", "M+15.995GLAM+15.995GGGGGASFDR"};
    
  for (auto input : inputs) {
      cout << "original input: " << input << endl;
      cout << "\tReplaceIWithL: " << Utility::Commons::ReplaceIWithL(input) << endl;
      cout << "\tRemovePTM_C: " << Utility::Commons::ReplacePTMOnCWithText(input, "(C)") << endl;
      cout << "\tReplacePTM_M: " << Utility::Commons::ReplacePTMOnMWithText(input, "(O)") << endl;
  }
  return 0;
}
