import os
import re
import sys
from copy import deepcopy

def run_save_iqtree(alignment_file, data_type, model, outname, true_tree, threads):
    
    ### Optimize tree under model ###
    if not os.path.exists(outname + "_inferredtree.iqtree"):
        print("Inferring tree", model)
        os.system("iqtree -nt " + str(threads) + " -m " + model + " -s " + alignment_file + " -st " + data_type)
        os.system("mv " + alignment_file + ".log " + outname + "_inferredtree.log")
        os.system("mv " + alignment_file + ".treefile " + outname + "_inferredtree.treefile")
        os.system("mv " + alignment_file + ".iqtree " + outname + "_inferredtree.iqtree")
        os.system("rm " + alignment_file + ".*")    
    
    ### Optimize parameters on true tree ###
    if not os.path.exists(outname + "_optimizedtruetree.iqtree"):
        print("Optimzing true tree", model)
        os.system("iqtree -nt " + str(threads) + " -m " + model + " -s " + alignment_file + " -st " + data_type + " -te " + true_tree)
        os.system("mv " + alignment_file + ".log " + outname + "_optimizedtruetree.log")
        os.system("mv " + alignment_file + ".treefile " + outname + "_optimizedtruetree.treefile")
        os.system("mv " + alignment_file + ".iqtree " + outname + "_optimizedtruetree.iqtree")
        os.system("rm " + alignment_file + ".*")     
      


def main():
    name                = sys.argv[1]
    tree                = sys.argv[2] 
    repl                = sys.argv[3]
    fitted_tree_path    = sys.argv[4]
    alignment_path      = sys.argv[5]
    model_file          = sys.argv[6]
    threads             = sys.argv[7]
   
   
    ##pandit
    if tree == "NA":
        rawname = name
        ###### !!!!!!!!! m1 tree here !!!!!!!!!
        true_tree_file = fitted_tree_path + name + "_m1_inferredtree.treefile"
    else: ## simulation
        rawname        = name + "_" + tree + "_rep" + repl + "_AA"
        true_tree_file =  "simulations/true_trees/" + tree + ".tree"
    alignment_file = alignment_path + rawname + ".fasta"
  
    with open(model_file, "r") as f: 
        all_models = f.readlines()
    use_models = {}
    for line in all_models:
        if line.startswith(name + "," + tree + "," + repl):
            line2 = line.split(",")
            model = line2[4].strip()
            q = line2[5].strip()
            use_models[q] = model
        
    # quartile models    
    for modelquant in range(1,6): 
        outname = fitted_tree_path + rawname + "_m" + str(modelquant)
        run_save_iqtree(alignment_file, "AA", use_models[str(modelquant)], outname, true_tree_file, threads)
    ## poisson
    outname = fitted_tree_path + rawname + "_poisson"
    run_save_iqtree(alignment_file, "AA", "Poisson", outname, true_tree_file, threads)
             

main()
