out=io_test
src=io_test.cc

read -p "Parse with FAST_PARSE? (yes|no): " ans

if [ "$ans" == "yes" ]; then
  FAST_PARSE="-D FAST_PARSE"
  out=$out.fast.parse
  echo 'Compiling with FAST_PARSE for mz, inten'
elif [ "$ans" == "no" ]; then
  FAST_PARSE=""
  out=$out.standard.parse
  echo 'Compiling with STANDARD_PARSE for mz, inten'
else 
  echo 'Invalid input. Exit'
  exit
fi

g++ -std=c++11 -O3 -ffast-math $FAST_PARSE ../class/*.cc ../utility/*.cc $src -o $out
echo "successfully compiled into exe:" $out
