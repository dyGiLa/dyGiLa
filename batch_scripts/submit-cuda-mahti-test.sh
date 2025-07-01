#!/bin/bash

#SBATCH --job-name dyGiLa
#SBATCH --time=00:15:00

#SBATCH --partition=gputest
#SBATCH --account=project_2014552

#SBATCH --nodes=1
#SBATCH --gres=gpu:a100:2

#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1

#SBATCH -e dyGiLa.e%j
#SBATCH -o dyGiLa.o%j

## Please remember to load the environment your application may need.
## And use the variable $LOCAL_SCRATCH in your batch job script
## to access the local fast storage on each node.
export LOCAL_SCRATCH=/scratch/project_2014552/test-dyGiLa-v0.0.2

module load openmpi/4.1.2-cuda
module unload netlib-scalapack/2.1.0 fftw/3.3.10-mpi

dstats="./stats"
dxmls="./rank_xmls"
dxdmf="./xdmf"
dinsitu="./insitu"
dpio="./pio"
dpioCurrent="./pio_Current"

if [ ! -d $dstats ] && [ ! -d $dxmls ] && [ ! -d $dxdmf ] && [ ! -d $dinsitu ] && [ ! -d $dpio ] && [ ! -d $dpioCurrent ]; then
    mkdir stats rank_xmls xdmf insitu pio pio_Current;
fi    

srun /projappl/project_2014552/dyGiLa-v0.0.2/build/dyGiLa -i sim_config_dyGiLa-Langevin-blob-quench-Hfield-30mT.txt
