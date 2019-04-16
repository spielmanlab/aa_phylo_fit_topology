library(ape)

tnames <- c("andersen", "dosreis", "greenalga", "opisthokonta", "prum", "salichos", "ruhfel", "rayfinned", "spiralia")

for (t in tnames)
{

    tree <- read.tree(paste0(t, ".tree"))
    tree2 <- multi2di(tree)
    write.tree(tree2, paste0(t, "_resolved.tree"))
}
