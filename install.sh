#!/bin/bash

echo ""
echo "Installing msSLASH and naive spectral searching method..."

CDIR=`pwd`
TMPDIR="src/msSLASH/bash_scripts"
cd $TMPDIR

echo 'Step 1: compiling msSLASH'
bash compile_search_msSLASH.sh
echo ''
echo 'Step 2: compiling naive spectral searching method'
bash compile_search_bruteforce.sh

cd $CDIR

if test -d bin; then
  read -n1 -p 'A local bin dir already exists, overwrite? (Y/N)' booleanYorN
  case $booleanYorN in
   y|Y) echo "" ; rm -f -r bin ; mkdir bin ;;
   n|N) echo "" ; echo "Only replacing binary files in local bin dir" ;;
   *) echo "" ; echo "Invalid Input "; exit 1 ;;
  esac
else 
  mkdir bin
fi

mv $TMPDIR/bruteforce bin/
mv $TMPDIR/msSLASH bin/

echo ""
echo "Executables are now installed under bin/"
echo ""

exit 0

