## Relationship between model fit and inference accuracy in protein phylogenetics

Please direct all questions or comments to Stephanie via `spielman <AT> rowan <DOT> edu`.

### Contents

+ `scripts/` contains all code for generating data, results, including SLURM submission scripts for HPC

+ `simulations/` contains all simulated alignments, DMS preferences, and phylogenies used for simulation. Simulations created with `scripts/simulate.py` and associated `scripts/submit_simulations.sbatch`.

+ `selected_models_<pandit/simulation>/` contains all results from model selection using `ModelFinder` with option `-m TESTONLY` to mimic other standard softwards. Contents created with `scripts/submit_model_selection_<pandit/simulation>.sbatch`

+ `processed_model_selection/` contains CSV's of results from model selection, including scores for all models as well as the specific "quantile" models. Contents created with with `scripts/process_model_selection.sh` (calls `scripts/selected_models_to_csv.py` and `scripts/parse_selected_models.R`)

+ `fitted_trees_ufb_<pandit/simulation>/` contains all inferred phylogenies with IQTREE and various associated logfiles. Contents created with `scripts/infer_trees.py` and associated `scripts/submit_tree_inference_<pandit/simulation>.sbatch`.


+ `topology_tests_<pandit/simulation>/` contains results from AU (and SH, unused in MS) tests as produced by `IQTREE`. Contents created with `scripts/run_topology_tests.py`.

+ `results_analysis/` contains a post-processing R script for stats+viz, relevant CSV files to post-process, and directory of figures used in MS produced by said script.
	+ `results_analysis/rf_fit_simulation.csv`, `results_analysis/ufb_splits_simulation.csv`, `results_analysis/rf_pandit.csv`, and `results_analysis/fit_pandit.csv` created with `scripts/calculate_rf_ufb.py`
	+ `results_analysis/topology_tests_<pandit/simulation.csv` created with `scripts/parse_topology_tests.py`
