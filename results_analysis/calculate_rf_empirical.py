import os
import sys
import dendropy
from dendropy.calculate import treecompare
import pprint

sim_path       = "simulations/"
truepath       = sim_path + "true_trees/"
inferencepath  = "fitted_trees_empirical/"
dms_list       = ["NP", "LAC", "Gal4", "HA", "HIV"]
treenames      = ["andersen", "dosreis", "greenalga", "opisthokonta", "prum", "ruhfel", "salichos", "rayfinned", "spiralia"]       
reps           = 20


fileinfo_order = ["name", "tree", "rep", "model", "optim"]
fitinfo_order  = ["logl", "k", "AIC", "AICc", "BIC"]

outfile = "inference_results_empirical.csv"
outstring = ",".join(fileinfo_order) + "," + ",".join(fitinfo_order) +",rf_true,treelength\n"

iqfiles_all = [x for x in os.listdir(inferencepath) if x.endswith(".iqtree")]

for dms in dms_list: 
    print(dms)
    for repindex in range(1,reps+1): ## files indexed from 1
        print(repindex)
        treespaces = {}
        truetrees = {}
        for treename in treenames:
            ts = dendropy.TaxonNamespace()
            treespaces[treename] = ts
            truetrees[treename] = dendropy.Tree.get_from_path(truepath + treename +"_resolved.tree", "newick", taxon_namespace = treespaces[treename], rooting='force-unrooted')
        
        for this_treename in treenames:
            print("  ",this_treename)
            prefix = dms + "_" + this_treename.replace(".tree","") + "_rep" + str(repindex) + "_AA"
            iqfiles = [x for x in iqfiles_all if x.startswith(prefix)]
            
            for iqfile in iqfiles:
                name = iqfile.replace(".iqtree","")
                treefile = name + ".treefile"

                fileinfo = {"name": None, "treerep": None, "tree": None, "model": None, "optim": None}
                components  = name.split("_")
                fileinfo["name"]  = components[0]
                fileinfo["tree"]  = components[1]
                fileinfo["rep"]   =  components[2].replace("rep","")
                fileinfo["model"]   = components[4]
                fileinfo["optim"]   = components[5]

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
            
with open(outfile, "w") as f:
    f.write(outstring.strip())
    
    
