library(cowplot)
library(tidyverse)
library(lme4)
library(lmerTest)
library(multcomp)
library(broom)
library(xtable)
library(colorspace)
library(ggridges)
library(ggforce)

theme_set(theme_classic() + theme(axis.line = element_line(colour = "grey10"),
                                  strip.background = element_rect(size=0.5)))
figure_directory <- "figures/"

########################## Factor levels and labeling ################################ 
name_levels <- c("1RII", "1R6M", "1IBS", "NP", "HA", "HIV")
name_labels_nsites <- c("1RII (195)", "1R6M (203)", "1IBS (291)", "NP (497)", "HA (564)", "HIV (661)")

model_levels <- c("m1", "m2", "m3", "m4", "m5", "poisson", "GTR20")
model_labels <- c("m1", "m2", "m3", "m4", "m5", "JC", "GTR")
m1to5_cols <- sequential_hcl(5, palette = "ylorrd")
poisson_col <- "grey40"
gtr20_col   <- "grey80" 
model_colors <- c(m1to5_cols, poisson_col, gtr20_col)
model_colors_nom1 <- c(m1to5_cols[2:5], poisson_col, gtr20_col)

model_levels_nom1 <- c("m2", "m3", "m4", "m5", "poisson", "GTR20")
model_labels_nom1 <- c("m2", "m3", "m4",  "m5", "JC", "GTR")

tree_levels <- c("ruhfel", "rayfinned", "dosreis", "prum", "andersen", "spiralia", "opisthokonta", "salichos")
tree_labels <- c("Green Plant", "Ray-finned fish", "Mammals", "Aves", "Lassa Virus", "Spiralia", "Opisthokonta", "Yeast")
tree_labels_abbr <- c("Plant", "Fish", "Mammals", "Aves", "Lassa", "Spiralia", "Opis.", "Yeast")
tree_labels_ntaxa <- c("Green Plant (360)", "Ray-finned fish (305)", "Mammals (274)", "Aves (200)", "Lassa Virus (179)", "Spiralia (103)", "Opisthokonta (70)", "Yeast (23)")

################################### Read in all data ##################################

### RF and BIC results
simulation_rf_fit  <- read_csv("rf_fit_simulation.csv",guess_max = 10000)
pandit_rf  <- read_csv("rf_pandit.csv",guess_max = 10000) 
pandit_fit <- read_csv("fit_pandit.csv",guess_max = 10000) 

### basic information about datasets
sim_info    <- read_csv("../simulations/simulation_trees_ntaxa.csv") 
pandit_info <- read_csv("../pandit_aa_alignments/pandit_info.csv")

### Topology test results
simulation_topology <- read_csv("topology_tests_simulation.csv")
pandit_topology     <- read_csv("topology_tests_pandit.csv")

## actual selected models
msel_simulation <- read_csv("../processed_model_selection/quantile_model_selection_simulation.csv")
all_sel         <- read_csv("../processed_model_selection/all_model_selection_simulation.csv")
msel_pandit     <- read_csv("../processed_model_selection/quantile_model_selection_pandit.csv")


## Bootstrap analysis
sim_ufb <- read_csv("ufb_splits_simulation.csv")

sim_ufb %>% 
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
        tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_ntaxa),
        name_levels  = factor(name, levels=name_levels), 
        rep = factor(rep)) -> ufb_fact

######### Normalize RF values #########
simulation_rf_fit %>% 
    left_join(sim_info) %>%
    mutate(max_rf = 2*ntaxa - 6) %>%  ### twice the number of internal edges
    mutate(rf_true_norm = rf/max_rf) -> simulation_rf_fit

pandit_rf %>% 
    left_join(pandit_info) %>%
    mutate(max_rf = 2*ntaxa - 6) %>%
    mutate(rf     = rf/max_rf) -> pandit_rf
pandit_fit %>% left_join(pandit_info) -> pandit_fit

##### entropy (for distinguishing simulations) ######
entropy <- read_csv("simulation_site_entropy.csv") %>%
                mutate(name_levels  = factor(name, levels=name_levels))



