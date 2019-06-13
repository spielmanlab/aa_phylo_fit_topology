library(tidyverse)

ufb <- read_csv("simulation_ufb_splits.csv")

ufb$model <- factor(ufb$model, levels=c("m1", "m2", "m3", "m4", "m5", "poisson", "GTR20"))
ufb$tree <- factor(ufb$tree)
ufb$rep <- factor(ufb$rep)
ufb$name <- factor(ufb$name)
ufb %>% 
    group_by(tree, model, rep) %>%
    tally() %>%
    rename(total_considered = n) -> ufb_total_nodes
     


############## False positive nodes ##############
ufb %>% 
    mutate(supported = boot >= 95) %>% 
    filter(in_true == FALSE, supported == TRUE) %>% 
    group_by(tree, rep, model) %>% 
    tally() %>% complete(tree, rep, model, fill = list(n = 0))%>%
    left_join(ufb_total_nodes) %>%
    mutate(percent = n / total_considered) %>%
    unique() %>%
    ggplot(aes(x = model, y = percent, color=model)) + 
        geom_jitter(width = 0.2, size=1) + 
        facet_wrap(~tree, scales="free_y", nrow=2) +
        geom_hline(yintercept = 0.05, color = "red") ### expected fp rate


############ False negative nodes ##############
ufb %>% 
    mutate(supported = boot >= 95) %>% 
    filter(in_true == TRUE, supported == FALSE) %>% 
    group_by(tree, rep, model) %>% 
    tally() %>% complete(tree, rep, model, fill = list(n = 0))%>%
    left_join(ufb_total_nodes) %>%
    mutate(percent = n / total_considered) %>%
    unique() %>%
    ggplot(aes(x = model, y = percent, color=model)) + 
        geom_jitter(width = 0.2, size=1) + 
        facet_wrap(~tree, scales="free_y", nrow=2)
        
        
#######################################################
######################## PR AUC #######################

library(PRROC)

final <- tibble(name = as.character(),
                tree = as.character(),
                model = as.character(),
                rep = as.numeric(),
                pr_auc = as.numeric())
                
for (n in unique(ufb$name)) {
    for (t in unique(ufb$tree)){
        print(t)
        for (m in unique(ufb$model)) {
            for (r in 1:20) {
                
                x <- ufb %>% 
                        filter(name == n, tree == t, model == m, rep == r) %>%
                        mutate(prob = boot / 100)
                fg <- x$prob[x$in_true == TRUE]
                bg <- x$prob[x$in_true == FALSE]
                pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
                this_pr_auc <- pr$auc.integral
                
                final <- bind_rows(final, 
                                    tibble(name = n, 
                                           tree = t, 
                                           model = m,
                                           rep = r, 
                                           pr_auc = this_pr_auc))
}}}}

                
                




























