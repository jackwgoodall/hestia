#!/bin/bash
#SBATCH --job-name=full_mod_transvir         
#SBATCH --mail-type=END,FAIL          
#SBATCH --mail-user=jack.goodall@lshtm.ac.uk 
#SBATCH --nodes=1                    
#SBATCH --ntasks=1                 
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=parallel_%j.log      
#SBATCH --cpus-per-task=52   

date;hostname;pwd

source ~/miniconda3/etc/profile.d/conda.sh

conda activate stan

echo "First full model of transvir... fingers crossed!"

Rscript hpc_first_transvir.R

date