#!/bin/bash


TYPE=$1
python3 processed_model_selection/selected_models_to_csv.py $TYPE ## -> all_model_selection<_empirical>.csv
Rscript processed_model_selection/parse_selected_models.R $TYPE  ## -> quantile_model_selection<_empirical>.csv