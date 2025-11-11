#!/bin/bash

dirarr=($(echo p-26.0-tauQ{50..3050..100}))
tauQarr=($(echo {50..3050..100}))

NumElem=${#tauQarr[@]}

for ((n=0;n<NumElem;n++)); do

    cd ${dirarr[$n]}

    pwd

    sbatch submit-ROCm.sh

    sleep 20s
    
    cd ..

done    

    

    



    


