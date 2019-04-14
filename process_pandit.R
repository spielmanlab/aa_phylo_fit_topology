library(cowplot)
library(tidyverse)
library(lme4)
library(multcomp)
library(broom)

theme_set(theme_classic())

pandit_tips <- read_csv("pandit_aa_alignments/info.csv") %>% 
                            mutate(max_rf = 2*nseq - 6) %>%
                            rename(name = pfam)
results <- read_csv("inference_results_pandit.csv") %>%
                        left_join(pandit_tips) %>%
                        mutate(rf_pandit_norm = rf_pandit / max_rf,
                               rf_q1_norm   = rf_q1 / max_rf)

model_levels <- c("pogofit", "q1", "q2", "q3", "q4", "q5", "poisson")
model_levels_sansq1 <- c("pogofit", "q2", "q3", "q4", "q5", "poisson")
model_labels <- c("Self-trained", "Selected M1", "Selected M2", "Selected M3", "Selected M4", "Selected M5", "Poisson")
model_labels_sansq1 <- c("Self-trained",  "Selected M2", "Selected M3", "Selected M4", "Selected M5", "Poisson")
ic_levels <- c("BIC", "AIC", "AICc")
################## pogo        m1        m2                
seven_colors <- c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4")
six_colors <- c("#d73027", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4")

results %>% 
  filter(optim == "inferredtree") %>%
  gather(ic, value, AIC, AICc, BIC) %>%
  dplyr::select(-k,-logl, -optim, -rf_pandit, -rf_pandit_norm, -rf_q1, -rf_q1_norm, -treelength) %>%
  group_by(ic, name) %>%
  mutate(ic.rank = as.integer(rank(value))) -> ic.ranks


### Fit rank with RANK on the x axis and fill by model. Much clearer than reverse.
ic.ranks %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
           ic_levels    = factor(ic, levels=ic_levels)) %>%
    ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black") + 
    facet_grid(~ic_levels) +
    scale_fill_manual(values = seven_colors, name = "Protein Model") +
    #geom_text(aes(label=model_levels),stat="count",position=position_stack(vjust=0.5), size=2.5)+
    xlab("Model rank") + ylab("Count")
    

## Are the ranks differences meaningful here?
# ic.ranks %>%
#     ungroup() %>%
#     group_by(name, ic) %>%
#     mutate(minIC = min(value), diffIC = abs(minIC - value)) %>%
#     ggplot(aes(x = ic.rank, y = diffIC, group = name)) + 
#         geom_point() + geom_line() + 
#         facet_wrap(~ic, scales = "free_y") + 
#         scale_x_continuous(breaks = 1:7, name = "Model rank") + 
#         panel_border() + background_grid() + 
#         ylab("Delta IC")



results %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
  ggplot(aes(x = model_levels, y = rf_pandit_norm,fill = model_levels)) + 
  geom_boxplot(outlier.shape = " ") + 
  geom_point(pch = 21, size=0.7, position = position_jitterdodge())+
  scale_fill_manual(values = seven_colors, name = "Protein Model") +
  ylab("RF Distance from PANDIT tree") + xlab("Protein Model") + 
  theme(legend.position = "none") +
  background_grid() -> p1
  
results %>%
  filter(optim == "inferredtree") %>% #, model != "q1") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
  ggplot(aes(x = model_levels, y = rf_q1_norm,fill = model_levels)) + 
  geom_boxplot(outlier.shape = " ") + 
  geom_point(pch = 21, size=0.7, position = position_jitterdodge())+
  scale_fill_manual(values = seven_colors, name = "Protein Model") +
  ylab("RF Distance from M1 tree") + xlab("Protein Model") + 
  theme(legend.position = "none") +
  background_grid() -> p2
  
plot_grid(p1, p2, nrow=2, labels="auto") -> pandit_rf_plots
  
  
  
sh <- read_csv("results_sh_pandit.csv")




