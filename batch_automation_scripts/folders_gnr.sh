#!/bin/bash

dirarr=($(echo p-26.0-tauQ{50..3050..100}))
tauQarr=($(echo {50..3050..100}))

NumElem=${#tauQarr[@]}

for ii in ${dirarr[@]}; do

    mkdir $ii && cd $ii

    pwd

    cp ../../test-dyGiLa-ROCm-standard-g/*.sh .
    cp ../../test-dyGiLa-ROCm-standard-g/*.txt .

    cd ..

done

for ((n=0;n<NumElem;n++)); do

    cd ${dirarr[$n]}

    pwd

    sed -i "s/tauQ1                50/tauQ1                ${tauQarr[${n}]}/g" sim_config_dyGiLa-Langevin.txt
    sed -i "s/tauQ2                50/tauQ2                ${tauQarr[${n}]}/g" sim_config_dyGiLa-Langevin.txt
    sed -i "s/seed                1/seed                 61/g" sim_config_dyGiLa-Langevin.txt
    # sed -i "s/IniT                 1.6687/IniT                 2.8536/g" sim_config_dyGiLa-Langevin*.txt    
    sed -i "s/Inip                26.0/Inip                6.0/g" sim_config_dyGiLa-Langevin*.txt;
    sed -i "s/InitH               0.0, 0.0, 0.0/InitH               0.0, 0.0, 150.0/g" sim_config_dyGiLa-Langevin.txt    
    # sed -i "s/do_gapA_clip         yes/do_gapA_clip         no/g" sim_config_dyGiLa-Langevin*.txt
    # sed -i "s/BCs1                    periodic/BCs1                    PairBreaking/g" sim_config_dyGiLa-Langevin*.txt
    # sed -i "s/BCs2                    periodic/BCs2                    PairBreaking/g" sim_config_dyGiLa-Langevin*.txt    

    # less sim_config_dyGiLa-Langevin*.txt
    
    cd ..

done    

    

    



    


