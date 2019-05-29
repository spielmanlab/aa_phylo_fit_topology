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

model_levels <- c("m1", "m2", "m3", "m4", "poisson", "m5")
model_labels <- c("M1", "M2", "M3", "M4", "Poisson", "M5")
model_levels_nom1 <- c("m2", "m3", "m4", "poisson", "m5")
model_labels_nom1 <- c("M2", "M3", "M4", "Poisson", "M5")

tree_levels <- c("ruhfel", "rayfinned", "dosreis", "prum", "andersen", "spiralia", "opisthokonta", "salichos")
tree_labels <- c("Green Plant", "Ray-finned fish", "Mammals", "Aves", "Lassa Virus", "Spiralia", "Opisthokonta", "Yeast")
tree_labels_twolines <- c("Green\nPlant", "Ray-finned\nfish", "Mammals", "Aves", "Lassa\nVirus", "Spiralia", "Opisthokonta", "Yeast")
tree_labels_ntaxa <- c("Green Plant (360)", "Ray-finned fish (305)", "Mammals (274)", "Aves (200)", "Lassa Virus (179)", "Spiralia (103)", "Opisthokonta (70)", "Yeast (23)")


################################### Read in all data ##################################

### RF and xIC results
simulation_rf_fit  <- read_csv("inference_results_simulation.csv",guess_max = 10000)
pandit_rf_fit     <- read_csv("inference_results_pandit.csv",guess_max = 10000) 

### basic information about datasets
sim_info <- read_csv("../simulations/simulation_trees_ntaxa.csv") 
pandit_info <- read_csv("../pandit_aa_alignments/info.csv")

### SH test results
simulation_sh <- read_csv("results_sh_simulation.csv")
pandit_sh <- read_csv("results_sh_pandit.csv")

## actual selected models
msel_simulation <- read_csv("../processed_model_selection/quantile_model_selection_simulation.csv")
all_sel <- read_csv("../processed_model_selection/all_model_selection_simulation.csv")
msel_pandit <- read_csv("../processed_model_selection/quantile_model_selection_pandit.csv")




######### Normalize RF values #########
simulation_rf_fit %>% 
    left_join(sim_info) %>%
    filter(tree != "greenalga") %>%
    mutate(max_rf = 2*ntaxa - 6) %>%
    mutate(rf_true_norm = rf_true/max_rf) -> simulation_rf_fit

pandit_rf_fit %>% 
    left_join(pandit_info) %>%
    mutate(max_rf = 2*ntaxa - 6) %>%
    mutate(rf_m1_norm     = rf_m1/max_rf) -> pandit_rf_fit


############################### Simulation figures ################################

#### simulation fit results
simulation_rf_fit %>%
  filter(optim == "inferredtree") %>% 
  dplyr::select(-AIC, -AICc, -k,-logl, -rf_true, -rf_true_norm, -treelength) %>%
  group_by(name, tree, rep) %>%
  mutate(ic.rank = as.integer(rank(BIC))) %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
           tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
           name_levels  = factor(name, levels=name_levels)) %>%
    ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black", size=.2) + 
    facet_grid(name_levels~tree_levels) +
    panel_border() + 
    scale_fill_brewer(palette = "RdYlBu", name = "Protein Model") +
    xlab("Model rank by BIC") + ylab("Count") + 
    theme(strip.text = element_text(size=8), 
          legend.position = "bottom", 
          legend.key.size = unit(4, "pt"), 
          legend.text = element_text(size=6), 
          legend.title = element_text(size=7), 
          axis.text = element_text(size=7), 
          axis.title = element_text(size=8)) +
    guides(fill = guide_legend(nrow=1)) -> model_fit_simulation_bars
l <- get_legend(model_fit_simulation_bars)
both <- plot_grid(model_fit_simulation_bars + theme(legend.position = "none"), l, nrow=2, rel_heights=c(1,0.05))
save_plot(paste0(figure_directory,"model_fit_simulation_bars.pdf"), both, base_width=7)
  
  
  
######## simulation RF Boxplots
simulation_rf_fit %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_twolines),
         name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
  ggplot(aes(x = model_levels, y = rf_true_norm, fill = model_levels)) + 
  geom_boxplot(outlier.size = 0.2, size=0.1) + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  facet_grid(tree_levels~name_levels, scales="free_y") +
  panel_border() +
  background_grid() +
  theme(axis.text.x = element_blank(),
       strip.text.y = element_text(size=8), 
       axis.ticks.x = element_blank(), 
       legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 1)) +
  xlab("Protein Models") + ylab("Normalized RF Distance") -> simulation_rf_boxplot
save_plot(paste0(figure_directory,"simulation_rf_boxplot_all.pdf"), simulation_rf_boxplot, base_height=8)


######## simulation RF Boxplots, NP only
simulation_rf_fit %>%
  filter(name == "NP", optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_ntaxa)) %>%
  ggplot(aes(x = model_levels, y = rf_true_norm, fill = model_levels)) + 
  geom_boxplot(outlier.size = 0.3, size=0.3) + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  facet_wrap(~tree_levels, scales="free_y", nrow = 2) +
  panel_border() +
  background_grid() +
  theme(axis.text.x = element_blank(),
       axis.text.y = element_text(size=7), 
       axis.ticks.x = element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1)) +
  xlab("Protein Models") + ylab("Normalized RF Distance") -> simulation_rf_boxplot_np
save_plot(paste0(figure_directory,"simulation_rf_boxplot_NP.pdf"), simulation_rf_boxplot_np, base_width = 10)





################################## PANDIT figures ###########################################


msel_pandit %>%
  filter(modelm == 1) %>%
  rowwise() %>%
  mutate(model_matrix = str_split(model_name, "\\+")[[1]][1]) %>%
  group_by(model_matrix) %>%
  tally() %>% 
  mutate(model_matrix_levels = fct_reorder(model_matrix, n)) %>%
  ggplot(aes(x = model_matrix_levels, y = n), fill = "grey60") + 
  geom_col() +
  geom_text(aes(x = model_matrix_levels, y = n+2, label = n), size=3)+
  xlab("M1 Model Matrix") + 
  ylab("Count") +
  theme(axis.text.x = element_text(size=8))-> m1_pandit_model_plot

  #### PANDIT fit results
  pandit_rf_fit %>%
    filter(optim == "inferredtree") %>%
    group_by(name) %>%
    mutate(ic.rank = as.integer(rank(BIC))) %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
    ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black", size=.2) + 
    panel_border() + 
    scale_fill_brewer(palette = "RdYlBu", name = "") +
    xlab("Model rank by BIC") + ylab("Count") + 
    theme(legend.position = "bottom") +
    guides( fill = guide_legend(nrow=1)) -> model_fit_pandit_bars
  
plot_grid(m1_pandit_model_plot, model_fit_pandit_bars, scale=0.93, nrow=1, labels="auto") -> pandit_model_grid

save_plot(paste0(figure_directory, "pandit_model_bars.pdf"), pandit_model_grid, base_width=10)
  



######## Pandit RF Boxplots
pandit_rf_fit %>%
  filter(model != "m1", optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels_nom1, labels = model_labels_nom1)) %>%
  ggplot(aes(x = model_levels, y = rf_m1_norm, fill = model_levels)) + 
  geom_boxplot() + 
  panel_border() + 
  scale_fill_manual(values = c("#FC8D59", "#FEE090", "#E0F3F8" ,"#91BFDB", "#4575B4"), name = "Protein Model", labels = model_labels_nom1) +
  #scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels_nom1) +
  xlab("Protein Models") + ylab("Normalized RF distance from M1 tree") +
  theme(legend.position = "none")-> rf_m1_pandit_plot     


  
#### PANDIT results
pandit_sh %>%
    gather(model, pvalue, m1:poisson) %>% 
    mutate(notsig = pvalue >= 0.01) %>%
    group_by(model, notsig) %>%
    tally() %>%
    mutate(perc_in_conf=n/200) %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) -> pandit_sh_summary

pandit_sh_summary  %>% filter(notsig == FALSE) %>% mutate(label = paste0(100*(perc_in_conf), "%")) -> pandit_sh_summary_false
pandit_sh_summary %>%
  mutate(notsig = factor(notsig, levels=c("TRUE", "FALSE"))) %>%
    ggplot(aes(x = model_levels, y = perc_in_conf, fill = notsig)) + 
        geom_col() +
        geom_text(data = pandit_sh_summary_false, aes(x = model_levels, y = perc_in_conf -.015, label = label), size=3) + 
        scale_fill_hue(l=60, name = "In M1 confidence set") +
        xlab("Models") + ylab("Percent of alignments") +
        theme(legend.position = "bottom") + 
        guides(fill = guide_legend(nrow=1, title.position = "left", hjust = 0.8))-> pandit_sh_barplot
        

pandit_sh %>% filter(m4 >=0.01, m5>=0.01, poisson <0.01) %>% nrow()  -> poisson_only   
pandit_sh %>% filter(m4 >=0.01, m5<0.01, poisson >=0.01) %>% nrow()  -> m5_only
pandit_sh %>% filter(m4 <0.01,  m5>=0.01, poisson >=0.01) %>% nrow() -> m4_only
pandit_sh %>% filter(m4 >=0.01, m5<0.01, poisson <0.01) %>% nrow()   -> m5_poisson       
pandit_sh %>% filter(m4 <0.01,  m5<0.01, poisson >=0.01) %>% nrow()  -> m4_m5
pandit_sh %>% filter(m4 <0.01,  m5>=0.01, poisson <0.01) %>% nrow()  -> m4_poisson
pandit_sh %>% filter(m4 <0.01,  m5<0.01, poisson <0.01) %>% nrow()   -> all_three        
venn_numbers <- tibble(type = c("m4_only", "m5_only", "poisson_only", "m4_m5", "m4_poisson", "m5_poisson", "all"), 
                      labels = c(m4_only, m5_only, poisson_only, m4_m5, m4_poisson, m5_poisson, all_three),                       x = c(0,     1.4,  -1.4,    0.75,   -0.75,   0,   0),   
                      y = c(1.5,  -0.7,   -0.7,    0.5,    0.5,    -1,   0))                                
df.venn <- data.frame(x = c(0, 0.866, -0.866),
                      y = c(1, -0.5, -0.5),
                      labels = c('M4', 'M5', 'Poisson'))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .4, size = 1, colour = 'grey60') +
  theme_void() + 
  coord_fixed(clip = "off") + 
theme(legend.position = 'none') +
  #scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
  #scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
  labs(fill = NULL) +
  annotate("text", x = venn_numbers$x, y = venn_numbers$y, label = venn_numbers$labels, size = 5) +
  annotate("text", x = 0, y = 2.65, label = "M4",fontface = "bold", size=4) +
  annotate("text", x = -2.6, y = -1.5,  label = "Poisson",fontface = "bold", size=4) +
  annotate("text", x = 2.45, y = -1.5, label = "M5",fontface = "bold", size=4) -> disagree_m1_venn




pandit_sh_barplot_nolegend <- pandit_sh_barplot + theme(legend.position = "none")
pandit_sh_barplot_legend<- get_legend(pandit_sh_barplot)
pandit_rf_sh_plot <- plot_grid(rf_m1_pandit_plot, pandit_sh_barplot_nolegend, disagree_m1_venn, scale=c(0.92, 0.92, 0.88), nrow=1, labels="auto")
pandit_rf_sh_plot_full <- plot_grid(pandit_rf_sh_plot, pandit_sh_barplot_legend, nrow=2, rel_widths=c(1, 0.2), rel_heights=c(1, 0.1))

save_plot(paste0(figure_directory,"pandit_rf_m1set.pdf"), pandit_rf_sh_plot_full, base_width=12)




######################### Model matrices and schematic ############################

#### Plot schematic for m1-m5
scheme_name <- "NP"
scheme_tree <- "opisthokonta"
scheme_rep <- 1
all_sel %>% 
  dplyr::select(name, tree, model, bic, repl) %>%
  filter(name == scheme_name, tree == scheme_tree, repl == scheme_rep) -> scheme_data

simulation_rf_fit %>% 
  filter(name == "NP", optim == "inferredtree", model == "poisson", rep == 1, tree == "opisthokonta") %>% 
  pull(BIC) -> poisson_bic


msel_simulation %>%
  filter(name == scheme_name, tree == scheme_tree, repl == scheme_rep) %>%
  rowwise() %>%
  mutate(modelm2 = paste0("M", modelm),
    model_levels_matrix = paste0(modelm2, " (", model_name, ")")) %>%
  bind_rows(
      tibble(name = scheme_name, tree = scheme_tree, repl = scheme_rep, modelm = 4.5, modelm2 = "Poisson", model_name = "Poisson", model_levels_matrix = "Poisson", bic = poisson_bic)) %>%
  mutate(model_levels= factor(modelm2, levels=c("M1", "M2", "M3", "M4", "Poisson", "M5")))   -> chunks


ggplot(scheme_data, aes(x = "", y = bic)) + 
    geom_violin(scale="width", fill="grey80", color = "grey70", alpha=0.6) + 
  geom_point(alpha = 0.3, size=2.5) +
  geom_point(data = chunks, shape=21, aes(x = 1, y = bic, fill = model_levels), size=2.5) + 
    geom_label(data = chunks, x = 1.02, hjust="outward", color = "black", aes(y = bic+250, fill = model_levels, label = model_levels_matrix), size=3, fontface = "bold", alpha=0.8) + 
 # scale_color_brewer(palette = "RdYlBu") +
  scale_fill_brewer(palette = "RdYlBu") +
  theme(legend.position = "none", axis.ticks.x = element_blank()) +
    xlab("") + ylab("BIC values across all tested models")  -> range_scheme_plot

ggsave(paste0(figure_directory, "bic_dist_scheme_qviolin.pdf"), range_scheme_plot, width = 6, height=5)
    

for (m in 1:5) {
    mname <- paste0("M", m)
    msel_simulation %>%
        filter(modelm == m) %>%
        rowwise() %>%
        mutate(model_matrix = str_split(model_name, "\\+")[[1]][1]) %>%
        group_by(name, tree, model_matrix) %>%
        tally() %>% 
        mutate(tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_twolines),
               name_levels  = factor(name, levels=name_levels),
               model_matrix_levels = fct_reorder(model_matrix, n)) %>%
        ggplot(aes(x = tree_levels, y = n, fill = model_matrix_levels)) + 
        geom_bar(stat="identity") +
        facet_wrap(~name_levels, nrow=2) + 
        xlab("Simulation tree") + 
        ylab("Count") + 
        labs(fill = paste(mname, "Model Matrix")) + 
        theme(axis.text.x = element_text(size=8)) -> selected_m_plot
    save_plot(paste0(figure_directory, "selected_", mname, "_simulation_barplot.pdf"), selected_m_plot, base_width=16, base_height = 6) 
}

      
    
x <- 1
plot_list <- c()
for (m in 1:5) {
    mname <- paste0("M", m)
    msel_pandit %>%
        filter(modelm == m) %>%
        rowwise() %>%
        mutate(model_matrix = str_split(model_name, "\\+")[[1]][1]) %>%
        group_by(model_matrix) %>%
        tally() %>% 
        mutate(model_matrix_levels = fct_reorder(model_matrix, n)) %>%
        ggplot(aes(x = model_matrix_levels, y = n, fill = model_matrix_levels)) + 
        geom_bar(stat="identity") +
        geom_text(aes(x = model_matrix_levels, y = n+2, label = n), size=3)+
        xlab("") + 
        ggtitle(paste(mname, "Model Matrix")) + 
        ylab("Count") + 
        #labs(fill = paste(mname, "Model Matrix")) + 
        theme(axis.text.x = element_text(size=8, angle=30, vjust=0.75),legend.position = "none") -> selected_m_plot
    if (m == 1) {
        selected_m_plot <- selected_m_plot + ggtitle("") + xlab(paste(mname, "Model Matrix")) 
        save_plot(paste0(figure_directory, "selected_", mname, "_pandit_barplot.pdf"), selected_m_plot, base_width=5, base_height = 3.25) 
    }
    if (m > 1) {
        plot_list[[x]] <- selected_m_plot
        x <- x+1
    } 
}
selected_grid <- plot_grid(plotlist = plot_list, nrow=2, labels = "auto")
save_plot(paste0(figure_directory, "selected_models_m2-5_pandit.pdf"), selected_grid, base_width=11, base_height=6) 





  
  


################################## SH tests ################################

#### simulation results
simulation_sh %>%
    gather(model, pvalue, m1:true) %>% 
    filter(pvalue < 0.01) %>%
    group_by(name, tree, model) %>%
    tally() %>% 
    arrange(model, name) %>%
    xtable()
# 1    NP opisthokonta   14    m5  0.009
# 2   HIV opisthokonta    4    m5  0.005
# 3   HIV opisthokonta   13    m5  0.006
# 4   HIV opisthokonta   15    m5  0.002
# 5   HIV       ruhfel   19    m5  0.001
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


# \begin{table}[ht]
# \centering
# \begin{tabular}{rlllr}
#   \hline
#  & name & tree & model & n \\ 
#   \hline
# 1 & HIV & opisthokonta & m5 &   3 \\ 
#   2 & HIV & ruhfel & m5 &   1 \\ 
#   3 & NP & opisthokonta & m5 &   1 \\ 
#   4 & Gal4 & andersen & true &   1 \\ 
#   5 & Gal4 & dosreis & true &   4 \\ 
#   6 & Gal4 & prum & true &   9 \\ 
#   7 & Gal4 & rayfinned & true &   8 \\ 
#   8 & Gal4 & ruhfel & true &  12 \\ 
#   9 & Gal4 & spiralia & true &   1 \\ 
#   10 & LAC & andersen & true &   1 \\ 
#    \hline
# \end{tabular}
# \end{table}


    
  



############################### Linear models ################################

######## RF for simulations
simulation_rf_fit %>% filter(optim == "inferredtree") -> simulation_rf_fit2
simulation_rf_fit2$model <- factor(simulation_rf_fit2$model, levels=c("m1", "m2", "m3", "m4", "m5", "poisson"))
lmer(rf_true_norm ~ model + (1|name) + (1|tree), data = simulation_rf_fit2) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
#                   Estimate Std. Error z value Pr(>|z|)  
# m2 - m1 == 0       0.004456   0.002654   1.679   0.5456    
# m3 - m1 == 0       0.007173   0.002654   2.703   0.0745 .  
# m4 - m1 == 0       0.008703   0.002654   3.280   0.0134 *  
# m5 - m1 == 0       0.026298   0.002654   9.911   <0.001 ***
# poisson - m1 == 0  0.003557   0.002654   1.340   0.7624    
# m3 - m2 == 0       0.002717   0.002654   1.024   0.9101    
# m4 - m2 == 0       0.004248   0.002654   1.601   0.5981    
# m5 - m2 == 0       0.021842   0.002654   8.231   <0.001 ***
# poisson - m2 == 0 -0.000899   0.002654  -0.339   0.9994    
# m4 - m3 == 0       0.001530   0.002654   0.577   0.9926    
# m5 - m3 == 0       0.019124   0.002654   7.207   <0.001 ***
# poisson - m3 == 0 -0.003616   0.002654  -1.363   0.7493    
# m5 - m4 == 0       0.017594   0.002654   6.631   <0.001 ***
# poisson - m4 == 0 -0.005147   0.002654  -1.940   0.3779    
# poisson - m5 == 0 -0.022741   0.002654  -8.570   <0.001 ***



######## RF against m1 tree for pandit
pandit_rf_fit %>% filter(optim == "inferredtree") -> pandit_rf_fit2
pandit_rf_fit2$model <- factor(pandit_rf_fit2$model, levels=c("m1", "m2", "m3", "m4", "m5", "poisson"))
lmer(rf_m1_norm ~ model + (1|name), data = pandit_rf_fit2) -> fit
glht(fit, linfct=mcp(model='Tukey')) %>% summary()
# Linear Hypotheses:
# m2 - m1 == 0       0.304420   0.008274  36.790  < 1e-04 ***
# m3 - m1 == 0       0.334983   0.008274  40.484  < 1e-04 ***
# m4 - m1 == 0       0.370579   0.008274  44.786  < 1e-04 ***
# m5 - m1 == 0       0.451461   0.008274  54.561  < 1e-04 ***
# poisson - m1 == 0  0.439134   0.008274  53.071  < 1e-04 ***
# m3 - m2 == 0       0.030563   0.008274   3.694 0.002962 ** 
# m4 - m2 == 0       0.066159   0.008274   7.996  < 1e-04 ***
# m5 - m2 == 0       0.147041   0.008274  17.771  < 1e-04 ***
# poisson - m2 == 0  0.134714   0.008274  16.281  < 1e-04 ***
# m4 - m3 == 0       0.035596   0.008274   4.302 0.000253 ***
# m5 - m3 == 0       0.116478   0.008274  14.077  < 1e-04 ***
# poisson - m3 == 0  0.104152   0.008274  12.587  < 1e-04 ***
# m5 - m4 == 0       0.080882   0.008274   9.775  < 1e-04 ***
# poisson - m4 == 0  0.068555   0.008274   8.285  < 1e-04 ***
# poisson - m5 == 0 -0.012327   0.008274  -1.490 0.670916    

    




######################### TREELENGTH FIGURES #################################

######## simulation treelength Boxplots
# simulation_rf_fit %>%
#   filter(optim == "inferredtree") %>%
#   mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
#          tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels_twolines),
#          name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
#   ggplot(aes(x = model_levels, y = treelength, fill = model_levels)) + 
#   geom_boxplot(outlier.size = 0.2, size=0.1) + 
#   scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
#   facet_grid(tree_levels~name_levels, scales="free_y") +
#   panel_border() +
#   background_grid() +
#   theme(axis.text.x = element_blank(),
#         strip.text.y = element_text(size=8), 
#         axis.ticks.x = element_blank(), 
#         legend.position = "bottom") +
#   guides(fill = guide_legend(nrow = 1)) +
#   xlab("Protein Models") + ylab("Inferred treelength") -> simulation_tl_boxplot
# save_plot(paste0(figure_directory,"simulation_tl_boxplot_all.pdf"), simulation_tl_boxplot, base_height=8)
# 
# 
# ######## simulation treelength Boxplots, NP only
# simulation_rf_fit %>%
#   filter(optim == "inferredtree") %>%
#   mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
#          tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels)) %>%
#   filter(name == "NP") %>%
#   ggplot(aes(x = model_levels, y = treelength, fill = model_levels)) + 
#   geom_boxplot(outlier.size = 0.3, size=0.1) + 
#   scale_fill_brewer(palette = "RdYlBu", name = "", labels = model_labels) +
#   facet_wrap(~tree_levels, scales="free_y", nrow = 2) +
#   panel_border() +
#   background_grid() +
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_text(size=7), 
#         axis.ticks.x = element_blank(), legend.position = "bottom") + 
#   guides(fill = guide_legend(nrow=1)) +
#   xlab("Protein Models") + ylab("Inferred treelength") -> simulation_tl_boxplot_np
# 
# 
# msel_simulation %>% 
#   filter(name == "NP", tree == "opisthokonta") %>% 
#   rowwise() %>% 
#   mutate(model_matrix = str_split(model, "\\+")[[1]][1],
#          model = paste0("r", modelr)) %>%
#   dplyr::select(-modelr) %>% 
#   rename(rep = repl) -> np_models
# 
# simulation_rf_fit %>% 
#   filter(name == "NP", optim == "inferredtree", tree == "opisthokonta") %>%
#   dplyr::select(name, tree, rep, treelength, model, BIC) %>%
#   left_join(np_models) %>%
#   mutate(model_matrix = ifelse(is.na(model_matrix), "Poisson", model_matrix)) %>%
#   mutate(model_levels = factor(model, levels = model_levels, labels = model_labels)) %>%
#   ggplot(aes(x = model_levels, y = treelength, group = rep)) + 
#   geom_point(aes(color = model_matrix), alpha = 0.9)  + geom_line(alpha=0.5, color="grey70") +
#   xlab("Protein Model") + ylab("Inferred treelength") + 
#   theme(legend.position = "bottom", legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.45, 'cm'), legend.key = element_rect(size=0.5), legend.text = element_text(size=6), legend.title= element_text(size=8) ) + 
#   panel_border() + 
#   guides(color = guide_legend(nrow=2)) + labs(color = "Model matrix") +
#   ggtitle("NP simulation replicates along the Opisthokonta tree") -> np_opis_tl_lineplot
# 
# bl_plots <- plot_grid(simulation_tl_boxplot_np, np_opis_tl_lineplot, nrow=1, labels="auto", scale=0.95)
# 
# save_plot(paste0(figure_directory,"simulation_tl_np_panels.pdf"), bl_plots, base_width = 12)
# 
# 
# 
# ############################### TL vs TL #################################
# 
# simulation_rf_fit %>% 
#   dplyr::select(name, tree, rep, model, optim, treelength) %>% 
#   spread(optim, treelength) %>% 
#   mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
#          tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
#          name_levels  = factor(name, levels=name_levels)) %>%
#   ggplot(aes(x = optimizedtruetree, y = inferredtree)) + 
#   geom_point(aes(color = tree_levels), alpha=0.8) + 
#   facet_grid(name_levels~model_levels) + 
#   geom_abline(color="red", alpha=0.7) +
#   coord_equal() + 
#   xlab("Treelength along optimized true tree") + 
#   ylab("Treelength along inferred model tree") +
#   labs(color = "Simulation tree") +
#   theme(axis.text = element_text(size=8))-> simulation_tl_inf_vs_optim
# ggsave(paste0(figure_directory, "simulation_tl_inf_vs_optim.pdf"), simulation_tl_inf_vs_optim, height = 5, width=7)
# 
# pandit_rf_fit %>% 
#   dplyr::select(name, model, optim, treelength) %>% 
#   spread(optim, treelength) %>% 
#   mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
#   ggplot(aes(x = optimizedtruetree, y = inferredtree)) + 
#   geom_point() + 
#   geom_abline(color="red") + 
#   coord_equal() + 
#   facet_wrap(~model_levels) + 
#   xlab("Treelength along optimized M1 tree") + 
#   ylab("Treelength along inferred model tree")  -> pandit_tl_inf_vs_optim
# ggsave(paste0(figure_directory, "pandit_tl_inf_vs_optim.pdf"), pandit_tl_inf_vs_optim,  height = 5, width=5)
# 
# 
# 
# 
# 
# #   
# #### PANDIT normalized treelength boxes
# pandit_rf_fit %>% 
#     filter(model == "m1", optim == "inferredtree") %>% 
#     dplyr::select(name, treelength) %>% 
#     rename(m1_treelength = treelength)  -> m1_tl
# pandit_rf_fit %>% 
#     filter(optim == "inferredtree") %>% 
#     dplyr::select(name, model, treelength) %>% 
#     left_join(m1_tl) %>% 
#     filter(model != "m1") %>%
#     mutate(tl_norm = treelength / m1_treelength,
#            model_levels = factor(model, levels=model_levels_nom1, labels = model_labels_nom1)) %>%
#     ggplot(aes(x = model_levels, y = tl_norm, fill = model_levels)) + 
#     geom_boxplot() + #color = "grey30", alpha=1, outlier.shape = " ") + 
#     panel_border() + 
#     #geom_jitter(alpha = 0.8, size=1) + 
#     geom_hline(yintercept=1, color = "#D73027", size=0.8) + 
#     scale_fill_manual(values = c("#FC8D59", "#FEE090", "#E0F3F8" ,"#91BFDB", "#4575B4"), name = "Protein Model", labels = model_labels_nom1) +
#     xlab("Protein Models") + ylab("Normalized treelength") +
#     theme(legend.position = "none") -> pandit_tl_norm





#### Box plot version  
# quant <- quantile(scheme_data$bic)
# scheme_data %>%
#   ggplot(aes(x = "", y = bic)) + 
#     geom_boxplot(fill = "firebrick3") +
#     xlab("") + ylab("BIC values across all tested models") +
#     theme(axis.ticks.x = element_blank(),
#           panel.border = element_blank(),
#           axis.line = element_line(color = 'black'), 
#           plot.margin = margin(0.5, 1, 0, 0.5, "cm")) +
#     geom_segment(x = 1.5, y = quant[[1]], xend = 1.02, yend = quant[[1]], color = "grey40", arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +
#     annotate(geom = "text", x = 1.55, y = min(scheme_data$bic), label = "M1") +
#     geom_segment(x = 1.5, y = quant[[2]], xend = 1.4, yend = quant[[2]], color = "grey40", arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +
#     annotate(geom = "text", x = 1.55, y = quant[[2]], label = "M2") +
#     geom_segment(x = 1.5, y = quant[[3]], xend = 1.4, yend = quant[[3]], color = "grey40", arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +
#     annotate(geom = "text", x = 1.55, y = quant[[3]], label = "M3") +
#     geom_segment(x = 1.5, y = quant[[4]], xend = 1.4, yend = quant[[4]], color = "grey40", arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +
#     annotate(geom = "text", x = 1.55, y = quant[[4]], label = "M4") +
#     geom_segment(x = 1.5, y = quant[[5]], xend = 1.02, yend = quant[[5]], color = "grey40", arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +
#     annotate(geom = "text", x = 1.55, y = quant[[5]], label = "M5")  -> scheme_plot
#    
# save_plot(paste0(figure_directory,"bic_dist_scheme.pdf"), scheme_plot) 