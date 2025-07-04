#!/bin/bash -l
#SBATCH --job-name=dyGILaG   # Job name
#SBATCH --output=dyGiLa-ROCm.o%j # Name of stdout output file
#SBATCH --error=dyGiLa-ROCm.e%j  # Name of stderr error file
#SBATCH --partition=dev-g  # partition name
#SBATCH --nodes=2               # Total number of nodes 
#SBATCH --ntasks-per-node=8     # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
#SBATCH --exclusive
#SBATCH --time=0-00:30:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000836  # Project for billing

# module --force purge && module --force unload LUMI
# module load LUMI/24.03 partition/G cpeCray/24.03 buildtools/24.03 rocm/6.0.3

# don't set SBATCH --exclusive if you use select_gpu

export MPICH_GPU_SUPPORT_ENABLED=1

dstats="./stats"
dxmls="./rank_xmls"
dxdmf="./xdmf"
dinsitu="./insitu"
dpio="./pio"
dpioCurrent="./pio_Current"

if [ ! -d $dstats ] && [ ! -d $dxmls ] && [ ! -d $dxdmf ] && [ ! -d $dinsitu ] && [ ! -d $dpio ] && [ ! -d $dpioCurrent ]; then
    mkdir stats rank_xmls xdmf insitu pio pio_Current;
fi    

srun /projappl/project_462000836/dyGiLa-IO/build/dyGiLa -i sim_config_dyGiLa-Langevin-blob-quench-Hfield-30mT.txt
#rm -rf ./select_gpu
