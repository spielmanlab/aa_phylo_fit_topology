library(ape)  


files <-c("anderson.tree", "opisthokonta.tree", "yeast.tree", "dosreis.tree", "prum.tree", "greenalga.tree", "ruhfel.tree", "greenplant.tree", "salichos.tree")


hi <- data.frame("tree" = as.character(), "ntaxa" = as.numeric())
for (tfile in files)
{
    t <- read.tree(tfile)
    t$edge.length <- t$edge.length * 3
    write.tree(t, paste0("times3/times3_", tfile))
    
    hi <- rbind(hi, data.frame("tree" = tfile, "ntaxa" = length(t$tip.label)))

}
write.csv(hi, "ntips.csv")


