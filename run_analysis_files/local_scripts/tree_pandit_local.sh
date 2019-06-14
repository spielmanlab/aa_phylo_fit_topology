#!/bin/bash

TOPPATH=/home/spielman/dms_modelselection
ALNPATH=${TOPPATH}/pandit_aa_alignments/
OUTPATH=${TOPPATH}/fitted_trees_pandit/

for LINE in {42..200}; do
    echo $LINE


    FILE=${ALNPATH}/names.txt
    INLINE=$(sed "${LINE}q;d" ${FILE})
    IFS=',' read -r NAME <<< "$INLINE"

    cd $TOPPATH
    python3 run_analysis_files/infer_trees.py $NAME NA NA $OUTPATH $ALNPATH processed_model_selection/quantile_model_selection_pandit.csv 4

done
