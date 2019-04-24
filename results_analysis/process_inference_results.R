library(cowplot)
library(tidyverse)
library(lme4)
library(lmerTest)
library(multcomp)
library(broom)
library(xtable)

theme_set(theme_bw())
figure_directory <- "figures/"

########################## Factor levels and labeling ################################ 
name_levels <- c("Gal4", "LAC", "NP", "HA", "HIV")
name_labels_nsites <- c("Gal4 (63)", "LAC (262)", "NP (497)", "HA (564)", "HIV (661)")

model_levels <- c("q1", "q2", "q3", "q4", "q5", "poisson")
model_labels <- c("M1", "M2", "M3", "M4", "M5", "Poisson")
model_levels_noq1 <- c("q2", "q3", "q4", "q5", "poisson")
model_labels_noq1 <- c("M2", "M3", "M4", "M5", "Poisson")

tree_levels <- c("ruhfel", "rayfinned", "dosreis", "prum", "andersen", "spiralia", "opisthokonta", "greenalga", "salichos")
tree_labels <- c("Green Plant", "Ray-finned fish", "Mammals", "Aves", "Lassa Virus", "Spiralia", "Opisthokonta", "Green Algae", "Yeast")
tree_labels_twolines <- c("Green\nPlant", "Ray-finned\nfish", "Mammals", "Aves", "Lassa\nVirus", "Spiralia", "Opisthokonta", "Green\nAlgae", "Yeast")
tree_labels_ntaxa <- c("Green Plant (360)", "Ray-finned fish (305)", "Mammals (274)", "Aves (200)", "Lassa Virus (179)", "Spiralia (103)", "Opisthokonta (70)", "Green Algae (23)", "Yeast (23)")


################################### Read in all data ##################################

### RF and xIC results
empirical_rf_fit  <- read_csv("inference_results_empirical.csv",guess_max = 10000)
pandit_rf_fit     <- read_csv("inference_results_pandit.csv",guess_max = 10000) 

### basic information about datasets
emp_info <- read_csv("../simulations/empirical_trees_ntaxa.csv") 
pandit_info <- read_csv("../pandit_aa_alignments/info.csv")

### SH test results
empirical_sh <- read_csv("results_sh_empirical.csv")
pandit_sh <- read_csv("results_sh_pandit.csv")



######### Normalize RF values #########
empirical_rf_fit %>% 
    left_join(emp_info) %>%
    mutate(max_rf = 2*ntaxa - 6) %>%
    mutate(rf_true_norm = rf_true/max_rf) -> empirical_rf_fit

pandit_rf_fit %>% 
    left_join(pandit_info) %>%
    mutate(max_rf = 2*ntaxa - 6) %>%
    mutate(rf_q1_norm     = rf_q1/max_rf) -> pandit_rf_fit




############################### Linear models ################################

######## Model One: RF for simulations

empirical_rf_fit %>% filter(optim == "inferredtree") -> empirical_inferred
empirical_inferred$model <- factor(empirical_inferred$model, levels=c("q1", "q2", "q3", "q4", "q5", "poisson"))
lmer(rf_true_norm ~ model + (1|name) + (1|tree), data = empirical_inferred) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
#                   Estimate Std. Error z value Pr(>|z|)  
# q2 - q1 == 0       0.004456   0.002654   1.679   0.5456    
# q3 - q1 == 0       0.007173   0.002654   2.703   0.0745 .  
# q4 - q1 == 0       0.008703   0.002654   3.280   0.0134 *  
# q5 - q1 == 0       0.026298   0.002654   9.911   <0.001 ***
# poisson - q1 == 0  0.003557   0.002654   1.340   0.7624    
# q3 - q2 == 0       0.002717   0.002654   1.024   0.9101    
# q4 - q2 == 0       0.004248   0.002654   1.601   0.5981    
# q5 - q2 == 0       0.021842   0.002654   8.231   <0.001 ***
# poisson - q2 == 0 -0.000899   0.002654  -0.339   0.9994    
# q4 - q3 == 0       0.001530   0.002654   0.577   0.9926    
# q5 - q3 == 0       0.019124   0.002654   7.207   <0.001 ***
# poisson - q3 == 0 -0.003616   0.002654  -1.363   0.7493    
# q5 - q4 == 0       0.017594   0.002654   6.631   <0.001 ***
# poisson - q4 == 0 -0.005147   0.002654  -1.940   0.3779    
# poisson - q5 == 0 -0.022741   0.002654  -8.570   <0.001 ***

######## Model Two: Treelength for simulations
empirical_rf_fit %>% filter(optim == "inferredtree") -> empirical_inferred
empirical_inferred$model <- factor(empirical_inferred$model, levels=c("q1", "q2", "q3", "q4", "q5", "poisson"))
lmer(treelength ~ model + (1|name) + (1|tree), data = empirical_inferred) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
#                   Estimate Std. Error z value Pr(>|z|)    
# q2 - q1 == 0        1.3669     0.2156   6.339   <1e-04 ***
# q3 - q1 == 0        1.0494     0.2156   4.866   <1e-04 ***
# q4 - q1 == 0        2.7004     0.2156  12.522   <1e-04 ***
# q5 - q1 == 0        5.0543     0.2156  23.438   <1e-04 ***
# poisson - q1 == 0  -1.7806     0.2156  -8.257   <1e-04 ***
# q3 - q2 == 0       -0.3175     0.2156  -1.472    0.682    
# q4 - q2 == 0        1.3335     0.2156   6.184   <1e-04 ***
# q5 - q2 == 0        3.6875     0.2156  17.100   <1e-04 ***
# poisson - q2 == 0  -3.1474     0.2156 -14.595   <1e-04 ***
# q4 - q3 == 0        1.6510     0.2156   7.656   <1e-04 ***
# q5 - q3 == 0        4.0050     0.2156  18.572   <1e-04 ***
# poisson - q3 == 0  -2.8300     0.2156 -13.123   <1e-04 ***
# q5 - q4 == 0        2.3540     0.2156  10.916   <1e-04 ***
# poisson - q4 == 0  -4.4809     0.2156 -20.779   <1e-04 ***
# poisson - q5 == 0  -6.8349     0.2156 -31.695   <1e-04 ***


######## Model Three: RF against q1 tree for pandit
pandit_rf_fit %>% filter(optim == "inferredtree") -> pandit_inferred
pandit_inferred$model <- factor(pandit_inferred$model, levels=c("q1", "q2", "q3", "q4", "q5", "poisson"))
lmer(rf_q1_norm ~ model + (1|name), data = pandit_inferred) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
# q2 - q1 == 0       0.304420   0.008274  36.790  < 1e-04 ***
# q3 - q1 == 0       0.334983   0.008274  40.484  < 1e-04 ***
# q4 - q1 == 0       0.370579   0.008274  44.786  < 1e-04 ***
# q5 - q1 == 0       0.451461   0.008274  54.561  < 1e-04 ***
# poisson - q1 == 0  0.439134   0.008274  53.071  < 1e-04 ***
# q3 - q2 == 0       0.030563   0.008274   3.694 0.002962 ** 
# q4 - q2 == 0       0.066159   0.008274   7.996  < 1e-04 ***
# q5 - q2 == 0       0.147041   0.008274  17.771  < 1e-04 ***
# poisson - q2 == 0  0.134714   0.008274  16.281  < 1e-04 ***
# q4 - q3 == 0       0.035596   0.008274   4.302 0.000253 ***
# q5 - q3 == 0       0.116478   0.008274  14.077  < 1e-04 ***
# poisson - q3 == 0  0.104152   0.008274  12.587  < 1e-04 ***
# q5 - q4 == 0       0.080882   0.008274   9.775  < 1e-04 ***
# poisson - q4 == 0  0.068555   0.008274   8.285  < 1e-04 ***
# poisson - q5 == 0 -0.012327   0.008274  -1.490 0.670916    


############################### Figures ################################

######## Empirical RF Boxplots
empirical_rf_fit %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
         name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
  ggplot(aes(x = model_levels, y = rf_true_norm, fill = model_levels)) + 
  geom_boxplot(outlier.shape = " ", size=0.1) + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  facet_grid(name_levels~tree_levels, scales="free_y") +
  panel_border() +
  background_grid() +
  theme(axis.text.x = element_blank(),
       axis.text.y = element_text(size=7), 
       axis.ticks.x = element_blank()) + 
  xlab("Protein Models") + ylab("Normalized RF Distance") -> empirical_rf_boxplot
save_plot(paste0(figure_directory,"empirical_rf_boxplot_all.pdf"), empirical_rf_boxplot, base_width = 10.5)


######## Empirical RF Boxplots, NP only
empirical_rf_fit %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_ntaxa)) %>%
  filter(name == "NP") %>%
  ggplot(aes(x = model_levels, y = rf_true_norm, fill = model_levels)) + 
  geom_boxplot(outlier.shape = " ", size=0.1) + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  facet_wrap(~tree_levels, scales="free_y", nrow = 3) +
  panel_border() +
  background_grid() +
  theme(axis.text.x = element_blank(),
       axis.text.y = element_text(size=7), 
       axis.ticks.x = element_blank()) + 
  xlab("Protein Models") + ylab("Normalized RF Distance") -> empirical_rf_boxplot_np
save_plot(paste0(figure_directory,"empirical_rf_boxplot_NP.pdf"), empirical_rf_boxplot_np, base_width = 8)




######## Empirical treelength Boxplots
empirical_rf_fit %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
         name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
  ggplot(aes(x = model_levels, y = treelength, fill = model_levels)) + 
  geom_boxplot(outlier.shape = " ", size=0.1) + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  facet_grid(name_levels~tree_levels, scales="free") +
  panel_border() +
  background_grid() +
  theme(axis.text.x = element_blank(),
       axis.text.y = element_text(size=7), 
       axis.ticks.x = element_blank()) + 
  xlab("Protein Models") + ylab("Inferred treelength") -> empirical_tl_boxplot
save_plot(paste0(figure_directory,"empirical_tl_boxplot_all.pdf"), empirical_tl_boxplot, base_width = 10.5)


######## Empirical treelength Boxplots, NP only
empirical_rf_fit %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_ntaxa)) %>%
  filter(name == "NP") %>%
  ggplot(aes(x = model_levels, y = treelength, fill = model_levels)) + 
  geom_boxplot(outlier.shape = " ", size=0.1) + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  facet_wrap(~tree_levels, scales="free_y", nrow = 3) +
  panel_border() +
  background_grid() +
  theme(axis.text.x = element_blank(),
       axis.text.y = element_text(size=7), 
       axis.ticks.x = element_blank()) + 
  xlab("Protein Models") + ylab("Inferred treelength") -> empirical_tl_boxplot_np
save_plot(paste0(figure_directory,"empirical_tl_boxplot_NP.pdf"), empirical_tl_boxplot_np, base_width = 8)



###########################################################################################
###########################################################################################


######## Pandit RF Boxplots
pandit_rf_fit %>%
  filter(optim == "inferredtree", model != "q1") %>%
  mutate(model_levels = factor(model, levels=model_levels_noq1, labels = model_labels_noq1)) %>%
  ggplot(aes(x = model_levels, y = rf_q1_norm, fill = model_levels)) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels_noq1) +
  xlab("Protein Models") + ylab("Normalized RF distance from M1 Tree") +
  theme(legend.position = "none")-> rf_q1_plot
        
save_plot(paste0(figure_directory,"pandit_rf_boxplots.pdf"), rf_q1_plot, base_width = 6)
  

sample_names <- base::sample( unique(pandit_rf_fit$name), 10 )
pandit_rf_fit %>%
  filter(optim == "inferredtree", model != "q1", name %in% sample_names) %>%
  mutate(model_levels = factor(model, levels=model_levels_noq1, labels = model_labels_noq1),
         name2 = paste0(name, " (", ntaxa, " taxa; ", nsites, " sites)")) %>%
  ggplot(aes(x = model_levels, y = rfw_q1)) + 
  geom_line(aes(group=name2), color = "grey40") + 
  geom_point(aes(color = model_levels), size=3) + 
  facet_wrap(~name2, nrow=2, scales="free_y") +
  xlab("Protein Models") + ylab("Weighted RF distance from M1 tree") +
  theme(legend.position = "none") -> rfw_q1_plot
save_plot(paste0(figure_directory,"pandit_wrf_lineplot_subset10.pdf"), rfw_q1_plot, base_width = 12)


pandit_rf_fit %>%
  filter(optim == "inferredtree", model != "q1") %>%
  mutate(model_levels = factor(model, levels=model_levels_noq1, labels = model_labels_noq1)) %>%
  ggplot(aes(x = model_levels, y = rfw_q1)) + 
  geom_line(aes(group=name), color = "grey40") + 
  geom_point(aes(color = model_levels), size=3) + 
  facet_wrap(~name, nrow=20, scales="free_y") +
  xlab("Protein Models") + ylab("Weighted RF distance from M1 tree") +
  theme(legend.position = "none") -> rfw_q1_plot_full
save_plot(paste0(figure_directory,"pandit_wrf_lineplot_allofthem.pdf"), rfw_q1_plot_full, base_height = 18)


pandit_rf_fit %>%
  filter(optim == "inferredtree", model != "q1") %>%
  dplyr::select(name, model, rfw_q1, BIC) %>%
  group_by(name) %>%
  mutate(wrf_rank = rank(rfw_q1),
         model_rank = rank(BIC)) -> pandit_wrf_ranks
lmerTest::lmer(wrf_rank ~ model_rank + (1|name), data = pandit_wrf_ranks) %>% summary()
#model_rank    0.62550    0.02470 998.00000   25.33   <2e-16 ***


        
        
######## Pandit tree length COV
pandit_rf_fit %>%
  filter(optim == "inferredtree") %>%
  group_by(name) %>%
  summarize(cov_treelength = sd(treelength)/mean(treelength)) %>%
  ggplot(aes(x = cov_treelength)) +
  geom_histogram(fill = "grey75", color = "black") +
  scale_x_continuous(breaks=seq(0.1, 0.8, 0.1)) + 
  scale_y_continuous(expand=c(0,0))+
  xlab("Coefficient of Variance of inferred tree lengths") +
  ylab("Count") -> cov_tl_histogram

  
save_plot(paste0(figure_directory,"cov_tl_histogram.pdf"), cov_tl_histogram)
  
  
###################### model fit ##########################################

#### Empirical fit results
empirical_rf_fit %>%
  filter(optim == "inferredtree") %>%
  gather(ic, value, AIC, AICc, BIC) %>%
  dplyr::select(-k,-logl, -optim, -rf_true, -rf_true_norm, -treelength) %>%
  group_by(ic, name, tree, rep) %>%
  mutate(ic.rank = as.integer(rank(value))) %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
           tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
           name_levels  = factor(name, levels=name_levels)) %>%
    filter(ic == "BIC") %>%
    ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black", size=.2) + 
    facet_grid(name_levels~tree_levels) +
    panel_border() + 
    scale_fill_brewer(palette = "RdYlBu", name = "Protein Model") +
    xlab("Model rank by BIC") + ylab("Count") + 
    theme(strip.text = element_text(size=9)) -> model_fit_empirical_bars
save_plot(paste0(figure_directory,"model_fit_empirical_bars.pdf"), model_fit_empirical_bars, base_width=10)
  

#### PANDIT fit results
pandit_rf_fit %>%
  filter(optim == "inferredtree") %>%
  gather(ic, value, AIC, AICc, BIC) %>%
  dplyr::select(-k,-logl, -optim, -rf_q1, -rf_q1_norm, -rfw_q1, -treelength) %>%
  group_by(ic, name) %>%
  mutate(ic.rank = as.integer(rank(value))) %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
    filter(ic == "BIC") %>%
    ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black", size=.2) + 
    panel_border() + 
    scale_fill_brewer(palette = "RdYlBu", name = "Protein Model") +
    xlab("Model rank by BIC") + ylab("Count") + 
    theme(strip.text = element_text(size=9)) -> model_fit_pandit_bars
save_plot(paste0(figure_directory,"model_fit_pandit_bars.pdf"), model_fit_pandit_bars, base_width=6)
  



################################## SH tests ################################

#### Empirical results
empirical_sh %>%
    gather(model, pvalue, q1:true) %>% 
    filter(pvalue <= 0.01) %>%
    group_by(name, tree, model) %>%
    tally() %>% 
    arrange(model, name) %>%
    xtable()
# 1    NP opisthokonta   14    q5  0.009
# 2   HIV opisthokonta    4    q5  0.005
# 3   HIV opisthokonta   13    q5  0.006
# 4   HIV opisthokonta   15    q5  0.002
# 5   HIV       ruhfel   19    q5  0.001
# 6   LAC     andersen    3  true  0.003
# 7  Gal4     andersen   17  true  0.010
# 8  Gal4      dosreis    6  true  0.010
# 9  Gal4      dosreis    7  true  0.005
# 10 Gal4      dosreis   13  true  0.009
# 11 Gal4      dosreis   16  true  0.005
# 12 Gal4         prum    2  true  0.005
# 13 Gal4         prum    3  true  0.002
# 14 Gal4         prum    5  true  0.005
# 15 Gal4         prum    8  true  0.006
# 16 Gal4         prum    9  true  0.010
# 17 Gal4         prum   13  true  0.005
# 18 Gal4         prum   15  true  0.009
# 19 Gal4         prum   17  true  0.008
# 20 Gal4         prum   20  true  0.009
# 21 Gal4       ruhfel    2  true  0.009
# 22 Gal4       ruhfel    3  true  0.000
# 23 Gal4       ruhfel    4  true  0.010
# 24 Gal4       ruhfel    5  true  0.004
# 25 Gal4       ruhfel    6  true  0.006
# 26 Gal4       ruhfel    9  true  0.004
# 27 Gal4       ruhfel   10  true  0.005
# 28 Gal4       ruhfel   11  true  0.010
# 29 Gal4       ruhfel   15  true  0.008
# 30 Gal4       ruhfel   18  true  0.003
# 31 Gal4       ruhfel   19  true  0.002
# 32 Gal4       ruhfel   20  true  0.000
# 33 Gal4    rayfinned    3  true  0.001
# 34 Gal4    rayfinned    4  true  0.004
# 35 Gal4    rayfinned    6  true  0.005
# 36 Gal4    rayfinned    8  true  0.002
# 37 Gal4    rayfinned    9  true  0.006
# 38 Gal4    rayfinned   14  true  0.003
# 39 Gal4    rayfinned   17  true  0.001
# 40 Gal4    rayfinned   19  true  0.003
# 41 Gal4     spiralia    9  true  0.008



#### PANDIT results
pandit_sh %>%
    gather(model, pvalue, q1:pandit) %>% 
    mutate(sig = pvalue <= 0.01) %>%
    group_by(name, sig) %>%
    tally() %>%
    filter(sig==TRUE) %>%
    left_join(pandit_info) %>%
    mutate(size = nsites*ntaxa) %>%
    ggplot(aes(x = size, y=n))+geom_point(alpha=0.5)
# 1 pandit     24
# 2 poisson    70
# 3 q4          7
# 4 q5         80
## 12% of the time differs pandit differs signidicantly
## 40% q5 differs


##################### Descriptive plot of where models came from ################


#### Plot schematic for q1-q5
msel <- read_csv("../processed_model_selection/quantile_model_selection_empirical.csv")
all_sel <- read_csv("../processed_model_selection/all_model_selection_empirical.csv")

scheme_name <- "HA"
scheme_tree <- "spiralia"
scheme_rep <- 1
all_sel %>% 
  dplyr::select(name, tree, model, bic, repl) %>%
  filter(name == scheme_name, tree == scheme_tree, repl == scheme_rep) %>%
  ggplot(aes(x = "", y = bic)) + 
    geom_boxplot(fill = "firebrick3") +
    xlab("") + ylab("BIC Distribution from ModelFinder") +
    theme(axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'))-> p
save_plot(paste0(figure_directory,"bic_dist_raw.pdf"), p) ## annotated in Graphic


msel %>%
    filter(modelq == 1) %>%
    rowwise() %>%
    mutate(model_matrix = str_split(model, "\\+")[[1]][1]) %>%
    group_by(name, tree, model_matrix) %>%
    tally() %>% 
    mutate(tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_twolines),
           name_levels  = factor(name, levels=name_levels),
           model_matrix_levels = factor(model_matrix, levels=c("JTT", "HIVb", "WAG"))) %>%
    ggplot(aes(x = tree_levels, y = n, fill = model_matrix_levels)) + 
    geom_bar(stat="identity") +
    facet_wrap(~name_levels) + 
    xlab("Simulation tree") + 
    ylab("Count") + 
    labs(fill = "M1 Model Matrix") + 
    theme(axis.text.x = element_text(size=7)) -> selected_m1_plot
save_plot(paste0(figure_directory,"selected_m1_empirical_barplot.pdf"), selected_m1_plot, base_width=16, base_height = 6) 
    
    

        
    
    
    
    
    
    
    
    
    
    
# empirical_rf_fit %>%
#   filter(optim == "inferredtree") %>%
#   gather(ic, value, AIC, AICc, BIC) %>%
#   dplyr::select(-k,-logl, -optim, -rf_true, -rf_true_norm, -treelength) %>%
#   group_by(ic, name, tree, rep) %>%
#   mutate(ic.rank = as.integer(rank(value))) -> ic.ranks
# 
# 
# ### Fit rank with RANK on the x axis and fill by model. Much clearer than reverse.
# ic.ranks %>%
#     mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
#            tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
#            name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
#     filter(ic == "BIC") %>%
#     ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
#     geom_bar(color="black", size=.2) + 
#     facet_grid(name_levels~tree_levels) +
#     panel_border() + 
#     scale_fill_brewer(palette = "RdYlBu", name = "Protein Model") +
#     xlab("Model rank by BIC") + ylab("Count") + 
#     theme(strip.text = element_text(size=9))
#     
# 
# ## Are the ranks differences meaningful here?
# ic.ranks %>%
#     ungroup() %>%
#     filter(ic == "BIC") %>%
#     group_by(name, tree, rep) %>%
#     mutate(minBIC = min(value), diffBIC = abs(minBIC - value)) %>%
#     ggplot(aes(x = ic.rank, y = diffBIC, group = rep)) + 
#         geom_point() + geom_line() + 
#         facet_wrap(name ~ tree, nrow=5, scales="free_y") + 
#         scale_x_continuous(breaks = 1:6, name = "Model rank") + 
#         panel_border() + background_grid() + 
#         ylab("Delta BIC")
    
    
    
