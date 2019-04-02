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


outstring = "name,rep,tree,model,optim,logl,AIC,BIC,AICc,k,logl,rf_true\n"

#NP_rtree64_bl3.0.tree_rep5_AA_q5_inferredtree.iqtree

fileinfo_order = ["name", "tree", "bl", "rep", "model", "optim"]
fitinfo_order  = ["logl", "k", "AIC", "AICc", "BIC"]




outstring = ",".join([x for x in fileinfo_order]) + ",rf_from_true," + ",".join([x for x in fitinfo_order]) + "\n"

iqfiles = [x for x in os.listdir(inferencepath) if x.endswith(".iqtree")]

for name in dms_list: 
    for this_treename in treenames:
        for repindex in range(1,6):
            prefix = name + "_" + this_treename + "_rep" + str(repindex) + "_AA"
            print(prefix)
            iqfiles = [x for x in os.listdir(inferencepath) if x.endswith(".iqtree") and x.startswith(prefix)]
        
            for iqfile in iqfiles:
                print(iqfile)
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
    

            # Log-likelihood of the tree: -6804.0639 (s.e. 261.9751)
            # Unconstrained log-likelihood (without tree): -1138.3914
            # Number of free parameters (#branches + #model parameters): 61
            # Akaike information criterion (AIC) score: 13730.1277
            # Corrected Akaike information criterion (AICc) score: 13779.5657
            # Bayesian information criterion (BIC) score: 13935.7367

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

                part = ",".join([fileinfo[x] for x in fileinfo_order]) + "," + rf + "," + ",".join([fitinfo[x] for x in fitinfo_order]) + "\n"
                print(part)
                outstring += part
                
    with open("inference_results.csv", "w") as f:
        f.write(outstring.strip())
    
    
