import os
import sys
import dendropy
from dendropy.calculate import treecompare
import pprint
truepath       = "simulations/true_trees/"


inferencepaths = {"rtree": "fitted_trees/", "empirical": "fitted_trees_empirical/"}
dms_list       = {"rtree": ["NP", "LAC", "Gal4", "HA", "HIV"], "empirical": ["HA"]}
treenames      = {"rtree": ["rtree100_bl0.03_rep1.tree", "rtree100_bl0.3_rep1.tree", "rtree100_bl0.75_rep1.tree", "rtree100_bl1.5_rep1.tree"], #,  "rtree100_bl3_rep1.tree"],
                  "empirical": ["anderson.tree", "dosreis.tree", "greenalga.tree", "greenplant.tree", "opisthokonta.tree", "prum.tree", "ruhfel.tree", "salichos.tree", "yeast.tree"]       
                 }
reps           = 20

def parse_components(type, fileinfo):
    if type  == "empirical":
        fileinfo["name"]  = components[0]
        fileinfo["tree"]  = components[1]
        fileinfo["rep"]   =  components[2].replace("rep","")
        fileinfo["model"]   = components[4]
        fileinfo["optim"]   = components[5]
    else:
        fileinfo["name"]    = components[0]
        fileinfo["tree"]    = components[2].split(".tree")[0].replace("bl", "")
        fileinfo["rep"]     = str(repindex)
        fileinfo["model"]   = components[6]
        fileinfo["optim"]   = components[7]
    return fileinfo


fileinfo_order = ["name", "tree", "rep", "model", "optim"]
fitinfo_order  = ["logl", "k", "AIC", "AICc", "BIC"]

outstring = "type," + ",".join(fileinfo_order) + "," + ",".join(fitinfo_order) +",rf_true,treelength\n"



for type in ["rtree", "empirical"]:
    iqfiles_all = [x for x in os.listdir(inferencepaths[type]) if x.endswith(".iqtree")]
    for dms in dms_list[type]: 
        for repindex in range(1,reps+1): ## files indexed from 1
            
            treespaces = {}
            truetrees = {}
            for treename in treenames[type]:
                ts = dendropy.TaxonNamespace()
                treespaces[treename] = ts
                truetrees[treename] = dendropy.Tree.get_from_path(truepath + treename , "newick", taxon_namespace = treespaces[treename], rooting='force-unrooted')
            
            for this_treename in treenames[type]:
            
                prefix = dms + "_" + this_treename.replace(".tree","") + "_rep" + str(repindex) + "_AA"
                iqfiles = [x for x in iqfiles_all if x.startswith(prefix)]
                
                for iqfile in iqfiles:
                    name = iqfile.replace(".iqtree","")
                    treefile = name + ".treefile"

                    fileinfo = {"name": None, "treerep": None, "tree": None, "model": None, "optim": None}
                    components  = name.split("_")
                    fileinfo = parse_components(type, fileinfo)
        

                    inftree = dendropy.Tree.get_from_path(
                                inferencepaths[type] + treefile,
                                "newick", 
                                taxon_namespace = treespaces[this_treename]
                               )
                    rf = str( treecompare.symmetric_difference( truetrees[this_treename], inftree) )
                    tl = str(inftree.length())
                
                    fitinfo = {"logl": None, "k": None, "AIC": None, "AICc": None, "BIC": None}
                    with open(inferencepaths[type] + iqfile, "r") as f:
                        iqlines = f.readlines()
                    for line in iqlines:
                        if line.startswith("Log-likelihood of the tree"):
                            fitinfo["logl"] = line.split(":")[1].split(" ")[1]
                        if line.startswith("Number of free parameters"):
                            fitinfo["k"] = line.split(" ")[-1].strip()
                        if line.startswith("Akaike information criterion"):
                            fitinfo["AIC"] = line.split(" ")[-1].strip()
                        if line.startswith("Corrected Akaike information criterion "):
                            fitinfo["AICc"] = line.split(" ")[-1].strip()
                        if line.startswith("Bayesian information criterion"):
                            fitinfo["BIC"] = line.split(" ")[-1].strip()
                    part = type + "," + ",".join([fileinfo[x] for x in fileinfo_order]) + "," + ",".join([fitinfo[x] for x in fitinfo_order]) + "," + rf + "," + tl  + "\n"
                    outstring += part
            
with open("inference_results.csv", "w") as f:
    f.write(outstring.strip())
    
    
