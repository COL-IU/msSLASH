#!/bin/bash


DIRECTORY=logs

if [ -d "$DIRECTORY" ]; then
    echo "clear $DIRECTORY"
    rm -rf $DIRECTORY
fi 
mkdir $DIRECTORY

#for cos_cut_off in `seq 0.4 0.1 0.9`; do
for cos_cut_off in `seq 0.5 0.1 0.5`; do
    echo 'working on cos_cut_off' $cos_cut_off
    for precision in `seq 0.1 0.1 1`; do
        echo '  precision: ' $precision
        #nohup time ./cluster-no-io $cos_cut_off $precision > $DIRECTORY/log_cos_${cos_cut_off}_pre_${precision} &
        #nohup time ./cluster-no-io-step20 $cos_cut_off $precision > $DIRECTORY/log_cos_${cos_cut_off}_pre_${precision} &
        nohup time ./cluster-no-io-step20-i200 $cos_cut_off $precision > $DIRECTORY/log_cos_${cos_cut_off}_pre_${precision} &
        #nohup time ./cluster-no-io-step20-i100 $cos_cut_off $precision > $DIRECTORY/log_cos_${cos_cut_off}_pre_${precision} &
        #nohup time ./cluster-no-io-step20-i50 $cos_cut_off $precision > $DIRECTORY/log_cos_${cos_cut_off}_pre_${precision} &
    done

done
