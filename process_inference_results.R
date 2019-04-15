library(cowplot)
library(tidyverse)
library(lme4)
library(multcomp)
library(broom)


dms <- tibble(name = c("Gal4", "LAC", "NP", "HA", "HIV"), nsites=c(63, 262, 497, 564, 661))

emp_tips <- read_csv("empirical_trees_ntaxa.csv") %>% mutate(max_rf = 2*ntaxa - 6)


results    <- read_csv("inference_results.csv",guess_max = 10000) 
sh_results <- read_csv("results_sh_empirical.csv")

name_levels <- c("Gal4", "LAC", "NP", "HA", "HIV")
name_labels_nsites <- c("Gal4 (63)", "LAC (262)", "NP (497)", "HA (564)", "HIV (661)")
model_levels <- c("q1", "q2", "q3", "q4", "q5", "poisson")
model_labels <- c("Selected M1", "Selected M2", "Selected M3", "Selected M4", "Selected M5", "Poisson")
tree_levels  <- c(0.03, 0.3, 0.75, 1.5) ## -> 0.1, 0.25,  0.5, 1.0
tree_labels  <- c("Low divergence", "Medium divergence", "High divergence", "Extra-high divergence") 
 
 
results %>% 
  filter(optim == "inferredtree") %>%
  gather(ic, value, AIC, AICc, BIC) %>%
  dplyr::select(-k,-logl, -optim, -rf_true, -rf_true_norm, -treelength) %>%
  group_by(ic, name, tree, rep) %>%
  mutate(ic.rank = as.integer(rank(value))) -> ic.ranks


### Fit rank with RANK on the x axis and fill by model. Much clearer than reverse.
ic.ranks %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
           tree_levels  = factor(tree, levels=tree_levels),
           name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
    filter(ic == "BIC") %>%
    ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black") + 
    facet_grid(name_levels~tree_levels) +
    panel_border() + 
    scale_fill_brewer(palette = "RdYlBu", name = "Protein Model") +
    geom_text(aes(label=model_levels),stat="count",position=position_stack(vjust=0.5), size=2.5)+
    xlab("Model rank by BIC")
    

### Are the ranks differences meaningful here?
ic.ranks %>%
    ungroup() %>%
    filter(ic == "BIC") %>%
    group_by(name, tree, rep) %>%
    mutate(minBIC = min(value), diffBIC = abs(minBIC - value)) %>%
    ggplot(aes(x = ic.rank, y = diffBIC, group = rep)) + 
        geom_point() + geom_line() + 
        facet_wrap(name ~ tree, nrow=5, scales="free_y") + 
        scale_x_continuous(breaks = 1:8, name = "Model rank") + 
        panel_border() + background_grid() + 
        ylab("Delta BIC")


rf_models <- tibble(tree = as.numeric(), term = as.character(), adj.p.value = as.numeric(), comparison = as.character(), estimate = as.numeric())
for (x in tree_levels){
    results %>% filter(tree == x, optim == "inferredtree") -> results2
    fit <- aov(rf_true ~ model + name, data= results2)
    TukeyHSD(fit) %>%
        tidy() %>%
        dplyr::select(adj.p.value, comparison,  estimate, term) %>%
        mutate(tree = x) -> x
    rf_models <- bind_rows(rf_models, x)
}
rf_models %>%
    filter(term == "model", adj.p.value <= 0.01) %>%
    arrange(tree, adj.p.value) %>%
    print.data.frame() 

results %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
         name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
  ggplot(aes(x = model_levels, y = rf_true_norm, fill = model_levels)) + 
  #geom_violin(scale="width")+
  geom_boxplot(outlier.shape = " ") + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  geom_jitter(width=0.3, height=0, size=0.75, alpha=0.7) +
  #facet_wrap(name_levels~tree_levels, nrow=5, scales="free_y") +
  facet_wrap(tree_levels~name_levels, nrow=4, scales="free_y") +
  panel_border() +
  background_grid() +
  theme(axis.text.x = element_text(size=10)) +    #element_text(size=9, angle=25, margin = margin(t = 8))) + 
  xlab("Protein Models") + ylab("Normalized RF Distance")

results %>%
  mutate(tree_levels  = factor(tree, levels=tree_levels),
         model_levels = factor(model, levels=model_levels)) %>%
  ggplot(aes(x = model, y = treelength)) + 
  geom_point(aes(color = name)) + background_grid() +
  facet_wrap(~tree_levels, nrow=1, scales="free_y") 

# results %>%
#   group_by(name, tree, model, rep) %>%
#   dplyr::select(-k, -AIC, -AICc, -BIC, -rf_true, -treelength) %>%
#   spread(optim, logl) %>%
#   mutate(true_minus_inf = optimizedtruetree - inferredtree,
#          model_levels = factor(model, levels=model_levels),
#          tree_levels  = factor(tree, levels=tree_levels),
#          name_levels  = factor(name, levels=name_levels)) %>%
#   filter(tree_levels == "0.75") %>%
#   ggplot(aes(x = true_minus_inf)) + 
#   geom_histogram() + 
#   scale_color_brewer(palette = "RdYlBu", name = "Protein Model") +
#   geom_vline(xintercept=0) + 
#   facet_wrap(name~model_levels, nrow=5, scales="free_y")


sh_results_rtree %>% 
    filter(type == "rtree") %>%
    dplyr::select(-type) %>%
    gather(model, pvalue, pogofit:true) %>% 
    mutate(sig = pvalue <= 0.01) %>% 
    group_by(tree, model) %>% 
    tally(sig) %>%
    ungroup() %>%
    ggplot(aes(x = model, y = n)) + 
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_text(aes(label=n, y = n+0.5)) + 
    facet_wrap(~tree, nrow=1) +
    panel_border()
    
sh_results_emp %>% 
    gather(model, pvalue, pogofit:true) %>% 
    mutate(sig = pvalue <= 0.01) %>% 
    group_by(tree, model) %>% 
    tally(sig) %>%
    ungroup() %>%
    ggplot(aes(x = model, y = n)) + 
    geom_bar(stat="identity", position = position_dodge()) + 
    facet_wrap(~tree, nrow=1) +
    panel_border()
    
  # 
sh_results_emp %>%
  dplyr::select(-type) %>%
  gather(model, pvalue, pogofit:true) %>%
#  mutate(model_levels = factor(model, levels=model_levels),
#         tree_levels  = factor(tree, levels=tree_levels, labels = tree_labels),
#         name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
#  filter(tree == "0.75") %>%
  ggplot(aes(x = pvalue, fill = name)) + 
  geom_histogram(color = "black") + 
  facet_grid(tree~model, scales="free")


#######################################################################################

results <- full_results %>% 
            filter(type == "empirical") %>%
            left_join(emp_tips, by = "tree") %>%
            group_by(tree) %>%
            mutate(rf_true_norm = rf_true / max_rf)
            
results %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels, labels = model_labels)) %>%
  ggplot(aes(x = model_levels, y = rf_true_norm, fill = model_levels)) + 
  geom_boxplot(outlier.shape = " ") + 
  scale_fill_brewer(palette = "RdYlBu", name = "Protein Model", labels = model_labels) +
  geom_jitter(width=0.3, height=0, size=0.75, alpha=0.7) +
  #facet_wrap(name_levels~tree_levels, nrow=5, scales="free_y") +
  facet_wrap(~tree, scales="free_y") +
  panel_border()

results %>%
  ggplot(aes(x = model, y = treelength, color=model, group=rep)) + 
  geom_point() + geom_line()+
  facet_grid(name~tree, scales="free_y")

rf_models <- tibble(tree = as.character(), term = as.character(), adj.p.value = as.numeric(), comparison = as.character(), estimate = as.numeric())
for (x in unique(results$tree)){
    results %>% filter(tree == x, optim == "inferredtree") -> results2
    fit <- aov(rf_true ~ model, data= results2)
    TukeyHSD(fit) %>%
        tidy() %>%
        dplyr::select(adj.p.value, comparison,  estimate, term) %>%
        mutate(tree = x) -> x
    rf_models <- bind_rows(rf_models, x)
}
rf_models %>%
    filter(term == "model", adj.p.value <= 0.01) %>%
    arrange(comparison, tree, adj.p.value) %>%
    print.data.frame() 
#            tree  term  adj.p.value      comparison estimate
# 1  opisthokonta model 3.918558e-03 poisson-hbstyle      4.4
# 2         yeast model 5.103848e-03 poisson-hbstyle      6.0
# 3  opisthokonta model 5.346448e-03 poisson-pogofit      4.3
# 4         yeast model 9.944940e-04 poisson-pogofit      6.7
# 5         yeast model 9.766320e-03      q4-hbstyle      5.7
# 6         yeast model 2.044820e-03      q4-pogofit      6.4
# 7       dosreis model 7.650857e-04      q5-hbstyle      7.4
# 8    greenplant model 5.636676e-04      q5-hbstyle     17.1
# 9  opisthokonta model 2.558466e-10      q5-hbstyle      8.4
# 10        yeast model 1.679588e-10      q5-hbstyle     11.8
# 11      dosreis model 3.007011e-04      q5-pogofit      7.8
# 12   greenplant model 2.036624e-03      q5-pogofit     15.8
# 13 opisthokonta model 4.155732e-10      q5-pogofit      8.3
# 14        yeast model 1.396494e-11      q5-pogofit     12.5
# 15        yeast model 7.895491e-03      q5-poisson      5.8
# 16 opisthokonta model 3.544106e-06           q5-q1      6.3
# 17        yeast model 4.238660e-06           q5-q1      8.7
# 18 opisthokonta model 8.072473e-06           q5-q2      6.1
# 19        yeast model 3.220690e-05           q5-q2      8.0
# 20 opisthokonta model 1.482904e-03           q5-q3      4.7
# 21        yeast model 3.642012e-04           q5-q3      7.1
# 22 opisthokonta model 7.518990e-04           q5-q4      4.9
# 23        yeast model 4.081584e-03           q5-q4      6.1
# == ->
#            tree  term  adj.p.value      comparison estimate
# 1       dosreis model 3.007011e-04      q5-pogofit      7.8
# 2       dosreis model 7.650857e-04      q5-hbstyle      7.4
# 3    greenplant model 5.636676e-04      q5-hbstyle     17.1
# 4    greenplant model 2.036624e-03      q5-pogofit     15.8
# 5  opisthokonta model 2.558466e-10      q5-hbstyle      8.4
# 6  opisthokonta model 4.155732e-10      q5-pogofit      8.3
# 7  opisthokonta model 3.544106e-06           q5-q1      6.3
# 8  opisthokonta model 8.072473e-06           q5-q2      6.1
# 9  opisthokonta model 7.518990e-04           q5-q4      4.9
# 10 opisthokonta model 1.482904e-03           q5-q3      4.7
# 11 opisthokonta model 3.918558e-03 poisson-hbstyle      4.4
# 12 opisthokonta model 5.346448e-03 poisson-pogofit      4.3
# 13        yeast model 1.396494e-11      q5-pogofit     12.5
# 14        yeast model 1.679588e-10      q5-hbstyle     11.8
# 15        yeast model 4.238660e-06           q5-q1      8.7
# 16        yeast model 3.220690e-05           q5-q2      8.0
# 17        yeast model 3.642012e-04           q5-q3      7.1
# 18        yeast model 9.944940e-04 poisson-pogofit      6.7
# 19        yeast model 2.044820e-03      q4-pogofit      6.4
# 20        yeast model 4.081584e-03           q5-q4      6.1
# 21        yeast model 5.103848e-03 poisson-hbstyle      6.0
# 22        yeast model 7.895491e-03      q5-poisson      5.8
# 23        yeast model 9.766320e-03      q4-hbstyle      5.7







#### Plot schematic for q1-q5
msel <- read_csv("quantile_model_selection.csv")

subres <- results %>% filter(optim == "inferredtree") %>% dplyr::select(name, tree, model, BIC, rep)
msel %>% 
    filter(name == "HA", tree == "rtree100_bl0.75_rep1", repl==1) %>%
    mutate(tree = "0.75") %>%
    rename(rep=repl, aamodel = model, model = modelq) %>%
    mutate(model = paste0("q", as.character(model))) %>%
    left_join(subres) -> qmodels_bic
    
results %>% 
    filter(name == "HA", rep == 1, tree == 0.75, optim == "inferredtree") %>%
    ggplot(aes(x = "", y = BIC)) + 
    geom_boxplot(fill = "dodgerblue2") + 
    xlab("") + ylab("Distribution of Alignment BIC scores from ModelFinder")  -> p
   #  
# p + geom_label_repel(data = model_bic, aes(x = 0.7, y = bic, label = model), hjust = 0)
#     arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
#     force = 10,
#     xlim  = x_limits
#   ) +
#   
# p <- p + annotate("text", x = 1.4, y = model_bic$bic[model_bic$comp == "diff_q1"], label =  model_bic$model[model_bic$comp == "diff_q1"])
# p + annotate("segment", x = 1.35, xend = 1.01, y = model_bic$bic[model_bic$comp == "diff_q1"], yend = model_bic$bic[model_bic$comp == "diff_q1"], arrow=arrow(type = "closed"), linejoin='mitre')
#     
#     
    
    
    
    
    
    
    
