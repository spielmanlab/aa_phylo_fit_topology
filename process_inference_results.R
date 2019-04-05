library(cowplot)
library(tidyverse)
library(lme4)
library(multcomp)
library(broom)


dms <- tibble(name = c("Gal4", "LAC", "NP", "HA", "HIV"), nsites=c(63, 262, 497, 564, 661))
max_rf <- 194 ## 2n - 6

results <- read_csv("inference_results.csv") %>% 
                select(-treerep) %>% 
                mutate(rf_true_norm = rf_true / max_rf)

name_levels <- c("Gal4", "LAC", "NP", "HA", "HIV")
name_labels_nsites <- c("Gal4 (63)", "LAC (262)", "NP (497)", "HA (564)", "HIV (661)")
model_levels <- c("pogofit", "hbstyle", "q1", "q2", "q3", "q4", "q5", "poisson")
model_labels <- c("Self-trained", "HB-style", "Selected Q1", "Selected Q2", "Selected Q3", "Selected Q4", "Selected Q5", "Poisson")
tree_levels  <- c(0.3, 1.5, 3) ## -> 0.1, 0.5, 1.0
tree_labels  <- c("Low divergence", "Medium divergence", "High divergence") ## -> 0.1, 0.5, 1.0
 
 
 
results %>% 
  filter(optim == "inferredtree") %>%
  gather(ic, value, AIC, AICc, BIC) %>%
  dplyr::select(-k,-logl, -optim, -rf_true, -rf_true_norm, -treelength) %>%
  group_by(ic, name, bl, rep) %>%
  mutate(ic.rank = as.integer(rank(value))) -> ic.ranks

# ### Fit rank with MODEL on the x axis. Less clear than rank on the X axis.
# ic.ranks %>%
#     mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
#            tree_levels  = factor(bl, levels=tree_levels,),
#            name_levels  = factor(name, levels=name_levels)) %>%
#     filter(ic == "BIC") %>%
#     ggplot(aes(x = model_levels, fill = factor(ic.rank))) + geom_bar(color="black") + 
#     facet_grid(name_levels~tree_levels) +
#     panel_border() + 
#     scale_fill_brewer(palette = "RdYlBu", name = "Model rank by BIC") +
#     geom_text(aes(label=ic.rank),stat="count",position=position_stack(vjust=0.5), size=3)+
#     xlab("Protein Models") -> barfits


### Fit rank with RANK on the x axis and fill by model. Much clearer than reverse.
ic.ranks %>%
    mutate(model_levels = factor(model, levels=model_levels, labels = model_labels),
           tree_levels  = factor(bl, levels=tree_levels,),
           name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
    filter(ic == "BIC") %>%
    ggplot(aes(x = factor(ic.rank), fill = model_levels)) + 
    geom_bar(color="black") + 
    facet_grid(name_levels~tree_levels) +
    panel_border() + 
    scale_fill_brewer(palette = "RdYlBu", name = "Protein Model") +
    geom_text(aes(label=model_levels),stat="count",position=position_stack(vjust=0.5), size=3)+
    xlab("Model rank by BIC")
    

### Are the ranks differences meaningful here?
ic.ranks %>%
    ungroup() %>%
    filter(ic == "BIC") %>%
    group_by(name, bl, rep) %>%
    mutate(minBIC = min(value), diffBIC = abs(minBIC - value)) %>%
    #mutate(ic.rank = factor(ic.rank, levels=c(1:8))) %>%
    #mutate(bic_diff = abs(value) - lag(abs(value), default = first(abs(value)))) %>%
    ggplot(aes(x = ic.rank, y = diffBIC, group = rep)) + 
        geom_point() + geom_line() + 
        facet_wrap(name ~ bl, nrow=5, scales="free_y") + 
        scale_x_continuous(breaks = 1:8, name = "Model rank") + 
        panel_border() + background_grid() + 
        ylab("Delta BIC")


rf_models <- tibble(bl = as.numeric(), term = as.character(), adj.p.value = as.numeric(), comparison = as.character(), estimate = as.numeric())
for (x in c(0.3, 1.5, 3)){
    results %>% filter(bl == x, optim == "inferredtree") -> results2
    fit <- aov(rf_true ~ model + name, data= results2)
    print(tidy(fit))
    TukeyHSD(fit) %>%
        tidy() %>%
        dplyr::select(adj.p.value, comparison,  estimate, term) %>%
        mutate(bl = bl) -> x
    rf_models <- bind_rows(rf_models, x)
}
rf_models %>%
    filter(term == "model", adj.p.value <= 0.01) %>%
    arrange(adj.p.value) %>%
    print.data.frame() ## all bl = 3


results %>%
  filter(optim == "inferredtree") %>%
  mutate(model_levels = factor(model, levels=model_levels),
         tree_levels  = factor(bl, levels=tree_levels, labels = tree_labels),
         name_levels  = factor(name, levels=name_levels, labels = name_labels_nsites)) %>%
  ggplot(aes(x = model_levels, y = rf_true_norm, fill = model_levels)) + 
  geom_violin(scale="width")+
  geom_jitter(width=0.3, height=0, size=1, alpha=0.5) +
  #facet_wrap(name_levels~tree_levels, nrow=5, scales="free_y") +
  facet_grid(name_levels~tree_levels, scales="free") +
  panel_border() +
  theme(axis.text.x = element_text(angle=30))

results %>%
  mutate(tree_levels  = factor(bl, levels=tree_levels)) %>%
  ggplot(aes(x = model, y = treelength)) + 
  geom_jitter(aes(color = name)) + 
  facet_wrap(~tree_levels, nrow=1,scales="free_y")

results %>%
  group_by(name, bl, model, rep) %>%
  select(-k, -AIC, -AICc, -BIC, -rf_true, -rf_true_norm, -treelength) %>%
  spread(optim, logl) %>%
  mutate(true_minus_inf = optimizedtruetree - inferredtree,
         model_levels = factor(model, levels=model_levels),
         tree_levels  = factor(bl, levels=tree_levels),
         name_levels  = factor(name, levels=name_levels)) %>%
  filter(tree_levels == "0.3") %>%
  ggplot(aes(x = true_minus_inf)) + 
  geom_histogram() + 
  geom_vline(xintercept=0) + 
  facet_wrap(name~model_levels, nrow=5, scales="free_y")

  
