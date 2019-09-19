src=search.cc
res=exe.search.omp.bruteforce

read -n1 -p 'Compile w/ FAST_PARSE?  (Y/N)' booleanYorN

case $booleanYorN in
  y|Y) echo "" ; echo "w/ FAST_PARSE"; FAST_PARSE=FAST_PARSE; res=$res.fast.parse ;;
  n|N) echo "" ; echo "w/o FAST_ARSE"; FAST_PARSE=None; res=$res.std.parse ;;
  *) echo "" ; echo "Invalid Input "; exit 1 ;;
esac

g++ -std=c++11 -O3 -ffast-math -D __BRUTEFORCE__ -D $FAST_PARSE ../utility/*.cc ../class/*.cc commons.cc $src -o $res -fopenmp

echo 'exe:' $res
