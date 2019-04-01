library(tidyverse)

dat <- read_csv("model_selection.csv")

selected_models <- tibble(name = as.character(),
                          tree = as.character(),
                          repl = as.integer(),
                          model = as.character(),
                          modelq = as.integer())



for (X in 1:5){
    dat %>% 
        group_by(name, tree, repl) %>% 
        mutate(thisbic = quantile(bic)[X],
               diffbic = abs(bic - thisbic)) %>%
        filter(diffbic == min(diffbic)) %>%
        ungroup() %>%
        dplyr::select(name, tree, repl, model) %>%
        mutate(modelq = X) -> tempdat
    selected_models <- bind_rows(selected_models, tempdat)
}
