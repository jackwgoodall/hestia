#!/bin/bash
#SBATCH --job-name=full_mod         # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jack.goodall@lshtm.ac.uk # Where to send mail
#SBATCH --nodes=1                     # Run all processes on a single node
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=parallel_%j.log      # Standard output and error log
#SBATCH --cpus-per-task=52   # 4 chains × 13 threads

date;hostname;pwd

source ~/miniconda3/etc/profile.d/conda.sh

conda activate stan

echo "First full model... fingers crossed!"

Rscript hpc_sim_setup.R

Rscript hpc_run_full_mod.R

date