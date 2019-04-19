library(cowplot)
library(tidyverse)
library(lme4)
library(multcomp)
library(broom)

theme_set(theme_bw())
figure_directory <- "figures/"

########################## Factor levels and labeling ################################ 
name_levels <- c("Gal4", "LAC", "NP", "HA", "HIV")
name_labels_nsites <- c("Gal4 (63)", "LAC (262)", "NP (497)", "HA (564)", "HIV (661)")
model_levels <- c("q1", "q2", "q3", "q4", "q5", "poisson")
model_labels <- c("M1", "M2", "M3", "M4", "M5", "Poisson")
tree_levels <- c("ruhfel", "rayfinned", "dosreis", "prum", "andersen", "spiralia", "opisthokonta", "greenalga", "salichos")
tree_labels <- c("Green Plant", "Ray-finned fish", "Mammals", "Aves", "Lassa Virus", "Spiralia", "Opisthokonta", "Green Algae", "Yeast")
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
    mutate(rf_pandit_norm = rf_pandit/max_rf,
           rf_q1_norm     = rf_q1/max_rf) -> pandit_rf_fit




############################### Linear models ################################

######## Model One: RF for simulations

empirical_rf_fit %>% filter(optim == "inferredtree") -> empirical_inferred
empirical_inferred$model <- factor(empirical_inferred$model, levels=c("q1", "q2", "q3", "q4", "q5", "poisson"))
lmer(rf_true ~ model + (1|name) + (1|tree), data = empirical_inferred) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
#                   Estimate Std. Error z value Pr(>|z|)  
# q2 - q1 == 0        1.0400     2.3360   0.445   0.9978  
# q3 - q1 == 0        1.6978     2.3360   0.727   0.9787  
# q4 - q1 == 0        2.1622     2.3360   0.926   0.9400  
# q5 - q1 == 0        6.6689     2.3360   2.855   0.0492 *
# poisson - q1 == 0   0.6689     2.3360   0.286   0.9997  
# q3 - q2 == 0        0.6578     2.3360   0.282   0.9998  
# q4 - q2 == 0        1.1222     2.3360   0.480   0.9968  
# q5 - q2 == 0        5.6289     2.3360   2.410   0.1528  
# poisson - q2 == 0  -0.3711     2.3360  -0.159   1.0000  
# q4 - q3 == 0        0.4644     2.3360   0.199   1.0000  
# q5 - q3 == 0        4.9711     2.3360   2.128   0.2726  
# poisson - q3 == 0  -1.0289     2.3360  -0.440   0.9979  
# q5 - q4 == 0        4.5067     2.3360   1.929   0.3841  
# poisson - q4 == 0  -1.4933     2.3360  -0.639   0.9880  
# poisson - q5 == 0  -6.0000     2.3360  -2.568   0.1052  


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



######## Model Three: RF against pandit tree for pandit
pandit_rf_fit %>% filter(optim == "inferredtree") -> pandit_inferred
pandit_inferred$model <- factor(pandit_inferred$model, levels=c("q1", "q2", "q3", "q4", "q5", "poisson"))
lmer(rf_pandit ~ model + (1|name), data = pandit_inferred) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
#
# q2 - q1 == 0         1.760      1.198   1.469   0.6840    
# q3 - q1 == 0        -0.080      1.198  -0.067   1.0000    
# q4 - q1 == 0         5.250      1.198   4.382   <0.001 ***
# q5 - q1 == 0        13.860      1.198  11.570   <0.001 ***
# poisson - q1 == 0   13.500      1.198  11.269   <0.001 ***
# q3 - q2 == 0        -1.840      1.198  -1.536   0.6409    
# q4 - q2 == 0         3.490      1.198   2.913   0.0418 *  
# q5 - q2 == 0        12.100      1.198  10.101   <0.001 ***
# poisson - q2 == 0   11.740      1.198   9.800   <0.001 ***
# q4 - q3 == 0         5.330      1.198   4.449   <0.001 ***
# q5 - q3 == 0        13.940      1.198  11.637   <0.001 ***
# poisson - q3 == 0   13.580      1.198  11.336   <0.001 ***
# q5 - q4 == 0         8.610      1.198   7.187   <0.001 ***
# poisson - q4 == 0    8.250      1.198   6.887   <0.001 ***
# poisson - q5 == 0   -0.360      1.198  -0.301   0.9997  


######## Model Four: RF against q1 tree for pandit
pandit_rf_fit %>% filter(optim == "inferredtree") -> pandit_inferred
pandit_inferred$model <- factor(pandit_inferred$model, levels=c("q1", "q2", "q3", "q4", "q5", "poisson"))
lmer(rf_q1 ~ model + (1|name), data = pandit_inferred) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
#
# q2 - q1 == 0         84.25       4.27  19.732  < 0.001 ***
# q3 - q1 == 0         91.79       4.27  21.498  < 0.001 ***
# q4 - q1 == 0        100.36       4.27  23.506  < 0.001 ***
# q5 - q1 == 0        117.53       4.27  27.527  < 0.001 ***
# poisson - q1 == 0   115.41       4.27  27.031  < 0.001 ***
# q3 - q2 == 0          7.54       4.27   1.766  0.48792    
# q4 - q2 == 0         16.11       4.27   3.773  0.00234 ** 
# q5 - q2 == 0         33.28       4.27   7.795  < 0.001 ***
# poisson - q2 == 0    31.16       4.27   7.298  < 0.001 ***
# q4 - q3 == 0          8.57       4.27   2.007  0.33805    
# q5 - q3 == 0         25.74       4.27   6.029  < 0.001 ***
# poisson - q3 == 0    23.62       4.27   5.532  < 0.001 ***
# q5 - q4 == 0         17.17       4.27   4.021  < 0.001 ***
# poisson - q4 == 0    15.05       4.27   3.525  0.00560 ** 
# poisson - q5 == 0    -2.12       4.27  -0.497  0.99631 


######## Model Five: treelength for pandit
pandit_rf_fit %>% filter(optim == "inferredtree") -> pandit_inferred
pandit_inferred$model <- factor(pandit_inferred$model, levels=c("q1", "q2", "q3", "q4", "q5", "poisson"))
lmer(treelength ~ model + (1|name) , data = pandit_inferred) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
#
#                   Estimate Std. Error z value Pr(>|z|)    
# q2 - q1 == 0         4.598      4.961   0.927    0.940    
# q3 - q1 == 0         3.009      4.961   0.607    0.991    
# q4 - q1 == 0        28.156      4.961   5.675   <1e-04 ***
# q5 - q1 == 0        72.187      4.961  14.551   <1e-04 ***
# poisson - q1 == 0  -27.975      4.961  -5.639   <1e-04 ***
# q3 - q2 == 0        -1.589      4.961  -0.320    1.000    
# q4 - q2 == 0        23.558      4.961   4.749   <1e-04 ***
# q5 - q2 == 0        67.590      4.961  13.624   <1e-04 ***
# poisson - q2 == 0  -32.573      4.961  -6.566   <1e-04 ***
# q4 - q3 == 0        25.147      4.961   5.069   <1e-04 ***
# q5 - q3 == 0        69.178      4.961  13.944   <1e-04 ***
# poisson - q3 == 0  -30.984      4.961  -6.245   <1e-04 ***
# q5 - q4 == 0        44.031      4.961   8.875   <1e-04 ***
# poisson - q4 == 0  -56.131      4.961 -11.314   <1e-04 ***
# poisson - q5 == 0 -100.163      4.961 -20.189   <1e-04 ***




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
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
  ggplot(aes(x = model_levels, y = rf_q1_norm, fill = model_levels)) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  xlab("Protein Models") + ylab("Normalized RF (M1)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=8))-> rf_q1_plot

pandit_rf_fit %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
  ggplot(aes(x = model_levels, y = rf_pandit_norm, fill = model_levels)) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  xlab("Protein Models") + ylab("Normalized RF (PANDIT)") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=8)) -> rf_pandit_plot
        
pandit_rf_plot <- plot_grid(rf_q1_plot, rf_pandit_plot, ncol=1, labels = "auto")
save_plot(paste0(figure_directory,"pandit_rf_boxplots.pdf"), pandit_rf_plot, base_height = 4.5)
  
######## Pandit RF Boxplots
pandit_rf_fit %>%
  filter(optim == "inferredtree") %>%
  group_by(name) %>%
  mutate(sd_treelength = sd(treelength), 
         range_treelength = range(treelength)[2] - range(treelength)[1]) %>%
  dplyr::select(name, sd_treelength, range_treelength, nsites, ntaxa) %>%
  unique() %>%
  mutate(size = nsites*ntaxa) -> pandit_tl_summary

# pandit_tl_summary %>%  
#   ggplot(aes(x = size, y = range_treelength)) +
#   geom_point() +
#   scale_x_log10() + scale_y_log10() + 
#   geom_smooth(method = "lm") + 
#   xlab("Dataset size") + 
#   ylab("Range of inferred treelengths")-> range_tl_plot

pandit_tl_summary %>%  
  ggplot(aes(x = size, y = sd_treelength)) +
  geom_point() +
  scale_x_log10() + scale_y_log10() + 
  geom_smooth(method = "lm") + 
  xlab("Dataset size") + 
  ylab("Std. dev. of inferred treelengths") -> sd_tl_plot
save_plot(paste0(figure_directory,"pandit_treelength_sd_plot.pdf"), sd_tl_plot)
  
# tl_plots <- plot_grid(sd_tl_plot, range_tl_plot, labels = "auto")


################################## SH tests ################################

#### Empirical results
empirical_sh %>%
    gather(model, pvalue, q1:true) %>% 
    filter(pvalue <= 0.01) %>%
    print.data.frame()
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
    
    
    
