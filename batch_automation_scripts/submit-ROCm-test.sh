#!/bin/bash -l
#SBATCH --job-name=dyGILaG   # Job name
#SBATCH --output=dyGiLa-ROCm.o%j # Name of stdout output file
#SBATCH --error=dyGiLa-ROCm.e%j  # Name of stderr error file
#SBATCH --partition=standard-g  # partition name
#SBATCH --nodes=1               # Total number of nodes 
#SBATCH --ntasks-per-node=8     # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
#SBATCH --time=2-00:00:00       # Run time (d-hh:mm:ss)
#SBATCH --account=project_462000836  # Project for billing

module --force purge && module --force unload LUMI
module load LUMI/24.03 partition/G cpeCray/24.03 buildtools/24.03 rocm/6.0.3

# don't set SBATCH --exlusive if you use select_gpu
cat << EOF > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"

export MPICH_GPU_SUPPORT_ENABLED=1

dstats="./stats"
dxmls="./rank_xmls"

if [ ! -d $dstats ] && [ ! -d $dxmls ]; then
    mkdir stats rank_xmls;
fi    

srun --cpu-bind=${CPU_BIND} ./select_gpu /projappl/project_462000836/dyGiLa-IO/build/dyGiLa -i sim_config_dyGiLa-Langevin-blob-quench-Hfield-30mT.txt
#srun /projappl/project_462000465/dyGiLa-ROCm/build/dyGiLa -i sim_config_dyGiLa-Langevin.txt
rm -rf ./select_gpu
