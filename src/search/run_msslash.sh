echo "Start time: " `date "+%Y%m%d-%Hh%Mm%Ss"`

while getopts t:d:i:c:e:p:m: option
do
case "${option}"
in
t) TARGET=${OPTARG};;
d) DECOY=${OPTARG};;
i) INTENSITY_SCALE_METHOD=${OPTARG};;
c) CHARGE=${OPTARG};;
e) EXP=${OPTARG};;
p) PSM=${OPTARG};;
m) MINZ=${OPTARG};;
esac
done

THREAD=10

echo -e "=== Target spectral library: " $TARGET
echo -e "=== Decoy spectral library: " $DECOY
echo -e "=== Charge: " $CHARGE
echo -e "=== Experimental spectral library: " $EXP
echo -e "=== Experimental spectral psm: " $PSM
echo -e "=== Min mz to consider: " $MINZ

if [ "$INTENSITY_SCALE_METHOD" == "raw" ]; then
  method=0
  echo "=== Use peak's raw intensity."
elif [ "$INTENSITY_SCALE_METHOD" == "log" ]; then
  method=1
  echo "=== LOG peak's intensity."
elif [ "$INTENSITY_SCALE_METHOD" == "sqrt" ]; then
  method=2
  echo "=== SQRT peak's intensity."
else
  echo "=== Unsupported method for peaks intensity. Exit.."
  exit 
fi
  
if [ $CHARGE == 2 ]; then
  precision=0.5
  echo "=== Fragment precision: " $precision
elif [ $CHARGE == 3 ]; then
  precision=0.33
  echo "=== Fragment precision: " $precision
else
  echo "=== Unsupported Charge. Exit.."
  exit 
fi

exp_prefix=${EXP/.mgf/}.peaks_${INTENSITY_SCALE_METHOD}
echo -e "=== Searching result prefix: "  $exp_prefix

echo -e "=== Searching with LSH\n"

for iter in 100;
do
  for h in 10;
  do
    out_file=./h${h}.i${iter}.precision${precision}.min_mz${MINZ}.threads${THREAD}.tsv
    echo $out_file
    /usr/bin/time -v ./exe.search.omp.lsh.fast.parse -n $h -i ${iter} -t ${THREAD} -l $TARGET -e $EXP -d $DECOY -o $out_file -r $method --precision $precision  --min_mz ${MINZ} -u
  done
done


echo -e "Searching with brute force\n"

out_file=$exp_prefix.bruteforce.precision${precision}.min_mz${MINZ}.threads${THREAD}.tsv
echo $out_file
/usr/bin/time -v ./exe.search.omp.bruteforce.fast.parse -t ${THREAD} -l $TARGET -e $EXP -d $DECOY -o $out_file -r $method --precision $precision  --min_mz ${MINZ} -u
