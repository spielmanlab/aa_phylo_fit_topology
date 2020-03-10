## Relative model fit does not predict topological accuracy in single-gene protein phylogenetics

Please direct all questions or comments to Stephanie via `spielman <AT> rowan <DOT> edu`.

### Contents

+ All directories and SLURM scripts are named as one of the following, representing the data they process:
  + `simulation`: MutSel simulations
  + `simulation_control`: WAG+I+G simulations
  + `pandit`: PANDIT analysis

+ `scripts/` contains all code for generating data, results, including SLURM (`*sbatch`) submission scripts for HPC

+ `simulations/` contains all simulated alignments, DMS preferences, and phylogenies used for simulation. Simulations were performed with `scripts/simulate.py` and `scripts/simulate_control.py` and associated `scripts/submit_simulations.sbatch` (for MutSel only, as those are slower than WAG which are quickly run locally).

+ `selected_models_*` contains all results from model selection using `ModelFinder` with option `-m TESTONLY`. Contents created with `scripts/submit_model_selection_*.sbatch`

+ `processed_model_selection/` contains CSV's of results from model selection, including scores for all models as well as the specific "quantile" models. Contents created with with `scripts/process_model_selection.sh` (calls `scripts/selected_models_to_csv.py` and `scripts/parse_selected_models.R`)

+ `fitted_trees_ufb_*` contains all inferred phylogenies with IQTREE and various associated logfiles. Contents created with `scripts/infer_trees.py` and associated `scripts/submit_tree_inference_*.sbatch`.

+ `topology_tests_*` contains results from AU tests as produced by `IQTREE`. Contents created with `scripts/run_topology_tests.py`.

+ `results_analysis/` contains a post-processing R script for stats+viz, relevant CSV files to post-process, and directory of figures used in MS produced by said script.
	+ `csv_files/` contains all processed CSVs created with one of the scripts in `scripts/`
	  + `all_models_pearson.csv` contain Pearson correlations among rate matrices (exchangabilities only!), created by `scripts/compare_all_models.py`
	  + `rf_fit_*.csv` contains robinson-foulds distances from true tree, and BIC, for all simulations. Contents created with `scripts/calculate_rf_ufb.py`.
	  + `rf_pandit.csv` contains all-to-all RF distances for PANDIT data. Contents created with `scripts/calculate_rf_ufb.py`.
	  + `topology_tests_*.csv` contains summarized results for topology tests. Contents created with `scripts/parse_topology_tests.py`
	+ `load.R` loads up data, libraries, and formats data for analysis in:
	  + `build_main_figures.R` creates all manuscript figures, saved in `main_figures/`
	  + `build_si_figures_tables.R` creates all SI manuscript figures and exports CSV tables, saved in `si_figures_tables/`
	  + `linear_model_and_misc_wrangling.Rmd` performs all linear models presented in the manuscript, and some miscellaneous data wrangling. 
	  
