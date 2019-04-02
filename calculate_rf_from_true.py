import os
import sys
import dendropy
from dendropy.calculate import treecompare
import pprint
inferencepath = "fitted_trees/"
truepath      = "simulations/true_trees/"
dms_list   = ["NP", "LAC", "Gal4", "HA", "HIV"]


treenames = [x for x in os.listdir(truepath) if x.endswith("tree")]
treespaces = {}
truetrees = {}
for treename in treenames:
    ts = dendropy.TaxonNamespace()
    treespaces[treename] = ts
    truetrees[treename] = dendropy.Tree.get_from_path(truepath + treename, "newick", taxon_namespace = treespaces[treename], rooting='force-unrooted')


fileinfo_order = ["name", "tree", "bl", "rep", "model", "optim"]
fitinfo_order  = ["logl", "k", "AIC", "AICc", "BIC"]

outstring = ",".join(fileinfo_order) + "," + ",".join(fitinfo_order) +",rf_true,treelength\n"


iqfiles_all = [x for x in os.listdir(inferencepath) if x.endswith(".iqtree")]

for dms in dms_list: 
    for this_treename in treenames:
        for repindex in range(1,6):
            prefix = dms + "_" + this_treename + "_rep" + str(repindex) + "_AA"
            iqfiles = [x for x in iqfiles_all if x.startswith(prefix)]

            for iqfile in iqfiles:
                name = iqfile.replace(".iqtree","")
                treefile = name + ".treefile"

                fileinfo = {"name": None, "tree": None, "bl": None, "model": None, "optim": None}
                #NP_rtree64_bl3.0.tree_rep5_AA_q5_inferredtree.iqtree
                components  = name.split("_")

                fileinfo["name"]    = components[0]
                fileinfo["tree"]    = components[1]
                fileinfo["bl"]      = components[2].split(".tree")[0].replace("bl", "")
                fileinfo["rep"]     = str(repindex)
                fileinfo["model"]       = components[5]
                fileinfo["optim"]       = components[6]

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
                part = ",".join([fileinfo[x] for x in fileinfo_order]) + "," + ",".join([fitinfo[x] for x in fitinfo_order]) + "," + rf + "," + tl  + "\n"
                outstring += part
            
with open("inference_results.csv", "w") as f:
    f.write(outstring.strip())
    
    
