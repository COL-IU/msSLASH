src=../search.cc
res=bruteforce
FAST_PARSE=FAST_PARSE
g++ -std=c++11 -O3 -ffast-math -D __BRUTEFORCE__ -D $FAST_PARSE ../../utility/*.cc ../../class/*.cc ../commons.cc $src -o $res -fopenmp
echo 'exe:' $res
