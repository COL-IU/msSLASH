#!/bin/bash

DIRECTORY=eval-res

if [ -d "$DIRECTORY" ]; then
    echo "clear $DIRECTORY"
    rm -rf $DIRECTORY
fi 
mkdir $DIRECTORY

ls out | grep cluster-final | while read -r line; do
    echo "evaluating $line"
    #nohup python eval.py out/$line > ${DIRECTORY}/res_${line} &
    python eval.py out/$line > ${DIRECTORY}/res_${line} 
done
