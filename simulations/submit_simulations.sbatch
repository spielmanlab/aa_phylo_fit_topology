#!/bin/bash
#SBATCH --job-name=simulate
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000
#SBATCH -p csm
#SBATCH --array=11-20

REP=${SLURM_ARRAY_TASK_ID} 

TOPDIR=/csm_data/spielman_lab/dms_modelselection/simulations
cd $TOPDIR
#mkdir -p alignments

python3 simulate_empirical.py $REP
#python3 simulate.py $REP
