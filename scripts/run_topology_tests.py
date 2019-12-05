import os
import sys
import re



model_order    = ["m1", "m2", "m3", "m4", "m5", "poisson", "GTR20"]

type                = sys.argv[1]
threads             = sys.argv[2]
treelist_name       = type + ".trees"
topology_test_path  = "../topology_tests_" + type + "/"
quantilefile        = "../processed_model_selection/quantile_model_selection_" + type + ".csv"
relll_reps          = "10000"

if type == "pandit":
    alignmentpath = "../pandit_aa_alignments/"
    inferencepath = "../fitted_trees_ufb_pandit/" 
elif type == "simulation":
    alignmentpath  = "../simulations/alignments/"
    inferencepath  = "../fitted_trees_ufb_simulation/"
    true_tree_path = "../simulations/true_trees/"
    sim_list       = ["LAC", "NP", "HA", "HIV"]
    treenames      = ["andersen", "dosreis", "opisthokonta", "prum", "ruhfel", "salichos", "rayfinned", "spiralia"]        
    reps           = 20


def call_test(threads, alnfile, modelname, treelist_name, relll_reps, m1tree, outfile):
    os.system("iqtree -quiet -nt " + threads + " -s " + alnfile + " -m " + modelname + " -z " + treelist_name + " -zb " + relll_reps + " -au -te " + m1tree + " -redo")
    os.system("mv " + alnfile + ".iqtree " + outfile)
    os.system("rm " + alnfile + ".*")

def determine_best_models(type, quantilefile):
    fitmodels = {}
    with open(quantilefile, "r") as f:
        qlines = f.readlines()
    
    if type == "pandit":
        for qline in qlines:
            name = qline.split(",")[0]
            model = qline.split(",")[2]
            quant = qline.split(",")[3].strip()
            if quant == "1":
                fitmodels[name] = model
    if type == "simulation":
        for qline in qlines:
            name = qline.split(",")[0]
            tree = qline.split(",")[1]
            rep = "rep" + qline.split(",")[2]
            prefix = "_".join([name, tree, rep, "AA"])
            model = qline.split(",")[4]
            quant = qline.split(",")[5].strip()
            if quant == "1":
                fitmodels[prefix] = model
    return(fitmodels)

def loop_over_tests(type):

    if type == "pandit":
        outstring = "name,m1,m2,m3,m4,m5,poisson,GTR20\n"
        for name in list(fitmodels.keys()): 
            print(name)     
            inferred_trees_raw = [x for x in os.listdir(inferencepath) if x.startswith(name) and x.endswith("inferredtree.treefile")]
            inferred_trees_ordered = [] 
            for model in model_order:
                for inftree in inferred_trees_raw:
                    if model in inftree:
                        with open(inferencepath + inftree, "r") as f:
                            ts = f.read().strip()
                        inferred_trees_ordered.append(ts)
                with open(treelist_name, "w") as f:
                    for ts in inferred_trees_ordered:
                        f.write(ts +"\n")
            m1tree = inferencepath + name + "_m1_inferredtree.treefile"
            alnfile = alignmentpath + name + ".fasta"
            
            ## Call the test
            outfile = topology_test_path + name + ".topology_tests"
            if not os.path.exists(outfile):
                call_test(threads, alnfile, fitmodels[name], treelist_name, relll_reps, m1tree, outfile)



    if type == "simulation":
        outstring  = "name,tree,repl,m1,m2,m3,m4,m5,poisson,GTR20,true\n"
        for name in sim_list:
            print(name)
            for tree in treenames:
                print("  ", tree)
                for rep in range(1,21):
                    print("  ", rep)
                    rawname = "_".join([name, tree, "rep" + str(rep), "AA"])
                    inferred_trees_raw = [x for x in os.listdir(inferencepath) if x.startswith(rawname) and x.endswith("inferredtree.treefile")]
                    inferred_trees_ordered = [] ### order based on replicate, model
                    for model in model_order:
                        for inftree in inferred_trees_raw:
                            if model in inftree:
                                with open(inferencepath + inftree, "r") as f:
                                    ts = f.read().strip()
                                inferred_trees_ordered.append(ts)
                    
                    with open(true_tree_path + tree + "_resolved.tree", "r") as f:
                        truetree = f.read().strip()
                    inferred_trees_ordered.append( truetree )
                    
                    with open(treelist_name, "w") as f:
                        for ts in inferred_trees_ordered:
                            f.write(ts +"\n")
                    
                    alnfile = alignmentpath + rawname + ".fasta"
                    m1tree = inferencepath + rawname + "_m1_inferredtree.treefile"
                    
                    ## Call the test
                    outfile = topology_test_path + rawname + ".topology_tests"
                    if not os.path.exists(outfile):
                        call_test(threads, alnfile, fitmodels[rawname], treelist_name, relll_reps, m1tree, outfile)
            
        
   
print("Determing all best-fitting models")
fitmodels = determine_best_models(type, quantilefile)
    
print("Run the tests")  
outstring = loop_over_tests(type)

