library(tidyverse)

dat <- read_csv("all_model_selection.csv")

selected_models <- tibble(name = as.character(),
                          tree = as.character(),
                          repl = as.integer(),
                          model = as.character(),
                          modelq = as.integer())




## gotta loop since sometimes multiple rows get returned for similarly fitting models
for (namex in unique(dat$name)){
    for(treex in unique(dat$tree)){
        for(replx in 1:20){
            dat %>% filter(name == namex, tree == treex, repl == replx) -> subdat
            for (q in 1:5) {
                subdat %>%
                    mutate(thisbic = quantile(bic)[q],
                           diffbic = abs(bic - thisbic)) %>%
                    filter(diffbic == min(diffbic)) %>%
                    ungroup() %>%
                    dplyr::select(name, tree, repl, model) %>%
                    mutate(modelq = q) -> tempdat
                    if (nrow(tempdat) >1 )
                    {
                        tempdat <- tempdat[1,]
                    }
                    selected_models <- bind_rows(selected_models, tempdat)
}}}}
write_csv(selected_models, "quantile_selected_models.csv")        
