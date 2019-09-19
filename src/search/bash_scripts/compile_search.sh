#!/bin/bash

#read -n1 -p 'Compile w/ FAST_PARSE?  (Y/N)' booleanYorN
#
#case $booleanYorN in
#  y|Y) echo "" ; echo "w/ FAST_PARSE"; op='y' ;;
#  n|N) echo "" ; echo "w/o FAST_ARSE"; op='n' ;;
#  *) echo "" ; echo "Invalid Input "; exit 1 ;;
#esac

echo 'Compiling searching with BF: ' 
bash /home/wang558/research/msms_analysis/src/search/bash_scripts/compile_search_bruteforce.sh $op

echo 'Compiling searching with LSH: ' 
bash /home/wang558/research/msms_analysis/src/search/bash_scripts/compile_search_lsh.sh $op
