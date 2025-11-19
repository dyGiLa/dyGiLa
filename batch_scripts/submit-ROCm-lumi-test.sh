#!/bin/bash -l
#SBATCH --job-name=dyGiLaG   # Job name
#SBATCH --output=dyGiLa-ROCm.o%j # Name of stdout output file
#SBATCH --error=dyGiLa-ROCm.e%j  # Name of stderr error file
#SBATCH --partition=dev-g # partition name
#SBATCH --nodes=2               # Total number of nodes 
#SBATCH --ntasks-per-node=8     # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
#SBATCH --exclusive
#SBATCH --time=0-00:30:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000960  # Project for billing

export MPICH_GPU_SUPPORT_ENABLED=1

dstats="./stats"
dxmls="./rank_xmls"
dxdmf="./xdmf"
dinsitu="./insitu"
dpio="./pio"
dpioCurrent="./pio_Current"

if [ ! -d $dstats ] && [ ! -d $dxmls ] && [ ! -d $dxdmf ] && [ ! -d $dinsitu ] && [ ! -d $dpio ] && [ ! -d $dpioCurrent ]; then
    mkdir stats rank_xmls xdmf insitu pio pio_Current;
    cd insitu && mkdir gapA-clip1 gapA-clip2 gapA-slice1 gapA-iso feDensity-clip_Camp Temeperature-clip Temperature-slice Temperature-iso pMarker-Slice pMarker-iso pMarker-fieldclip pMarker-isoVolume-Bphase pMarker-fieldclip-Aphase;
    cd ..
fi    

# srun --cpu-bind=${CPU_BIND} ./select_gpu /projappl/project_462000960/dyGiLa-develop/build/dyGiLa -i sim_config_dyGiLa-Langevin-blob-quench-Hfield-30mT.txt
srun /projappl/project_462000960/dyGiLa-develop-lite/build/dyGiLa -i sim_config_dyGiLa-Langevin-blob-quench-Hfield-30mT.txt

rm -rf ./select_gpu
