#!/bin/bash
#SBATCH --job-name=WT
#SBATCH --account=def-panch
#SBATCH --gpus=h100:1
#SBATCH --mem=20g
#SBATCH -t 168:00:00                 # time (days-hours:minutes:seconds)


module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.6 amber-pmemd/24.3

export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1

bash ./WT_run1/pmemd_cuda_job.sh &
bash ./WT_run2/pmemd_cuda_job.sh &
bash ./WT_run3/pmemd_cuda_job.sh &

wait
