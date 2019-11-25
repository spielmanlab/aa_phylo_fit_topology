#!/bin/bash
#SBATCH --job-name=simulate
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000
#SBATCH -p csm
#SBATCH --array=1-20

REP=${SLURM_ARRAY_TASK_ID} 

TOPDIR=/csm_data/spielman_lab/aa_phylo_fit_topology/scripts
cd $TOPDIR
#mkdir -p alignments

#python3 simulate_empirical.py $REP
for NAME in LAC HA NP HIV; do
  python3 simulate.py $NAME $REP
done
