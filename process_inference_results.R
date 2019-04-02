library(cowplot)
library(tidyverse)

dms <- tibble(name = c("Gal4", "LAC", "NP", "HA", "HIV"), nsites=c(63, 262, 497, 564, 661))
results <- read_csv("inference_results.csv") %>% unite(tree, tree, bl)

name_levels <- c("Gal4", "LAC", "NP", "HA", "HIV")
model_levels <- c("pogofit", "hbstyle", "poisson", "q1", "q2", "q3", "q4", "q5")
tree_levels  <- c("rtree64_0.3", "btree64_0.3", "rtree64_3", "btree64_3")

results %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels),
         tree_levels  = factor(tree, levels=tree_levels),
         name_levels  = factor(name, levels=name_levels)) %>%
  filter(name_levels == "HIV") %>%
  ggplot(aes(x = model_levels, y = rf_true, fill = model_levels)) + 
  geom_violin() + 
  facet_wrap(~tree_levels)

results %>%
  mutate(tree_levels  = factor(tree, levels=tree_levels)) %>%
  ggplot(aes(x = model, y = treelength)) + 
  geom_jitter(aes(color = name)) + 
  facet_wrap(~tree_levels, nrow=1,scales="free_y")

results %>%
  group_by(name, tree, model, rep) %>%
  select(-k, -AIC, -AICc, -BIC, -rf_true, -treelength) %>%
  spread(optim, logl) %>%
  mutate(true_minus_inf = optimizedtruetree - inferredtree,
         model_levels = factor(model, levels=model_levels),
         tree_levels  = factor(tree, levels=tree_levels),
         name_levels  = factor(name, levels=name_levels)) %>%
  filter(tree_levels == "rtree64_3") %>%
  ggplot(aes(x = true_minus_inf)) + 
  geom_histogram() + 
  geom_vline(xintercept=0) + 
  facet_wrap(name~model_levels, nrow=5, scales="free_y")

  

#   
#     spread(optim, logl) %>%
#  mutate(true_minus_optim = truetree - optimizedtree ) -> aa.optimvstrue
#   aa.optimvstrue %>%
#   filter(tree == "baltree128") %>%
#   mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
#   ggplot(aes(x = true_minus_optim, fill=model_levels)) +
#   geom_histogram(bins = 10, color="black") + 
#   geom_vline(xintercept =0,color="black", size=1) + 
#   facet_wrap(~model_levels, scales="free_x", nrow=1) +
#   scale_fill_hue(name = "Models", l=50) + 
#   xlab("True Tree LogL - Optimized Tree LogL") + 
#   ylab("Count") +
#   theme(legend.position = "none") -> logldiff_hist
# save_plot("logl_hist_n128.pdf", logldiff_hist, base_width = 10, base_height=2.75)
# 

