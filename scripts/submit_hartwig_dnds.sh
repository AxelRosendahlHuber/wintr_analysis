#!/bin/bash
#SBATCH -c 2
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --array=1-27%15
#SBATCH --job-name=HMF_wintr
#SBATCH --output="output/slurm-%j.out"
#SBATCH --partition=spot_cpu
#SBATCH --requeue

cd ~/wintr/hartwig_tests
source activate R 
Rscript ~/wintr/hartwig_tests/hartwig_dnds.R /home/arosendahl/wintr/hartwig_tests/hartwig_muts/ ${SLURM_ARRAY_TASK_ID} HMF
