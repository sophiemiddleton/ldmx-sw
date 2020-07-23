#!/bin/bash
source scripts/ldmx-env.sh
COUNTER=0
for filename in signalData/*; do
    ldmx fire MultiStage.py $filename
    mv output.root SignalOutput_$COUNTER.root
    let COUNTER+=COUNTER+1
done
