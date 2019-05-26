#!/bin/bash
#SBATCH --job-name=simulate
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID} 

TOPDIR=/csm_data/spielman_lab/dms_modelselection/simulations
cd $TOPDIR
#mkdir -p alignments

#python3 simulate_empirical.py $REP
for NAME in HA NP LAC Gal4 HIV; do
  python3 simulatepy $NAME $REP
done
