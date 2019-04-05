library(ape)  


files <-c("anderson.tree", "opisthokonta.tree", "yeast.tree", "dosreis.tree", "prum.tree", "greenalga.tree", "ruhfel.tree", "greenplant.tree", "salichos.tree")


for (tfile in files)
{

    t <- read.tree(tfile)
    t$edge.length <- t$edge.length * 3
    write.tree(t, paste0("times3_", tfile))
}