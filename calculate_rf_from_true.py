import os
import sys
import dendropy
from dendropy.calculate import treecompare
import pprint
inferencepaths = {"rtree": "fitted_trees/", "empirical": "fitted_trees_empirical/"}
truepath      = "simulations/true_trees/"
dms_list   = {"rtree": ["NP", "LAC", "Gal4", "HA", "HIV"], "empirical": "HA"}
reps = = {"rtree": 21, "empirical": 11}

treenames = {"rtree": ["rtree100_bl0.3_rep1.tree", "tree100_bl1.5_rep1.tree", "rtree100_bl0.75_rep1.tree", "rtree100_bl3_rep1.tree"],
             "empirical": ["anderson.tree", "dosreis.tree", "greenalga.tree", "greenplant.tree", "opisthokonta.tree", "prum.tree", "ruhfel.tree", "salichos.tree", "yeast.tree"]       
            }


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

outstring = "type,".join(fileinfo_order) + "," + ",".join(fitinfo_order) +",rf_true,treelength\n"


iqfiles_all = [x for x in os.listdir(inferencepath) if x.endswith(".iqtree")]

for type in ["rtree", "empirical"]:
    treespaces = {}
    truetrees = {}
    for treename in treenames[type]:
        ts = dendropy.TaxonNamespace()
        treespaces[treename] = ts
        truetrees[treename] = dendropy.Tree.get_from_path(truepath + treename + ".tree", "newick", taxon_namespace = treespaces[treename], rooting='force-unrooted')
        for this_treename in treenames:
            for dms in dms_list[type]: 
                for repindex in range(1,21):
                    prefix = dms + "_" + this_treename + "_rep" + str(repindex) + "_AA"
                    iqfiles = [x for x in iqfiles_all if x.startswith(prefix)]

                    for iqfile in iqfiles:
                        name = iqfile.replace(".iqtree","")
                        treefile = name + ".treefile"

                        fileinfo = {"name": None, "treerep": None, "tree": None, "model": None, "optim": None}
                        components  = name.split("_")
                        fileinfo = parse_components(type, fileinfo)
        

                        inftree = dendropy.Tree.get_from_path(
                                    inferencepath + treefile,
                                    "newick", 
                                    taxon_namespace = treespaces[this_treename]
                                   )
                        rf = str( treecompare.symmetric_difference( truetrees[this_treename], inftree) )
                        tl = str(inftree.length())
                
                        fitinfo = {"logl": None, "k": None, "AIC": None, "AICc": None, "BIC": None}
                        with open(inferencepath + iqfile, "r") as f:
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
    
    
