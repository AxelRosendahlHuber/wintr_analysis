#!/bin/bash
#SBATCH -c 2
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --array=1-27%15
#SBATCH --job-name=HMF_wintr
#SBATCH --output="output/slurm-%j.out"
#SBATCH --partition=spot_cpu
#SBATCH --requeue

cd /data/bbg/projects/mut_risk/wintr_analysis
source activate R

DATAFOLDER="/data/bbg/projects/mut_risk/raw_data/mutation_data/PCAWG_unfiltered/"
DATATYPE="PCAWG"

Rscript scripts/cohort_dnds.R -f ${DATAFOLDER} -i ${SLURM_ARRAY_TASK_ID} -n ${DATATYPE}
