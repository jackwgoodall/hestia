#!/bin/bash
#SBATCH --job-name=simulation         # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jack.goodall@lshtm.ac.uk # Where to send mail
#SBATCH --nodes=1                     # Run all processes on a single node
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --mem=32gb                     # Total memory limit
#SBATCH --time=02:00:00               # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log      # Standard output and error log
date;hostname;pwd

source ~/miniconda3/etc/profile.d/conda.sh

conda activate stan

echo "Simulation Set up"

Rscript hpc_sim_setup.R

date