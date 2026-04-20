#!/bin/bash
#SBATCH --job-name=simulation         # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jack.goodall@lshtm.ac.uk # Where to send mail
#SBATCH --nodes=1                     # Run all processes on a single node
#SBATCH --ntasks=4                   # Run a single task
#SBATCH --nodes=1                     # Keeps the job running on a single node.
#SBATCH --mem=1gb                     # Total memory limit
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log      # Standard output and error log
date;hostname;pwd

conda 

conda activate stan

echo "Simulation Set up"

Rscript sim_setup.R

date