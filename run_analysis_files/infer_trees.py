import os
import re
import sys
from copy import deepcopy

true_tree_path    = "simulations/true_trees/"
hb_path           = "hb_models/"


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
    name      = sys.argv[1]
    tree      = sys.argv[2] 
    repl      = sys.argv[3]
    fitted_tree_path   = sys.argv[4]
    alignment_path = sys.argv[5]
    pogomodel_path = sys.argv[6]
    model_file = sys.argv[7]
    threads   = sys.argv[8]
   
    
    if tree == "NA":
        rawname = name
        ###### !!!!!!!!! q1 tree here !!!!!!!!!
        true_tree_file = fitted_tree_path + name + "_q1_inferredtree.treefile"
    else:
        rawname        = name + "_" + tree + "_rep" + repl + "_AA"
        true_tree_file = true_tree_path + tree + ".tree"
    alignment_file = alignment_path + rawname + ".fasta"
  
    with open(model_file, "r") as f: 
        all_models = f.readlines()
    print(all_models)   
    use_models = {}
    # LAC,btree64_bl3.0.tree,5,HIVw,5
    print (name + "," + tree + "," + repl)
    for line in all_models:
        print(line)
        if line.startswith(name + "," + tree + "," + repl):
            line2 = line.split(",")
            model = line2[3].strip()
            q = line2[4].strip()
            use_models[q] = model
    print(use_models)
        
    # quartile models    
    for modelquant in range(1,6):
        modelquant = str(int_modelquant)
        outname = fitted_tree_path + rawname + "_q" + str(modelquant)
        run_save_iqtree(alignment_file, "AA", use_models[str(modelquant)], outname, true_tree_file, threads)
    ## poisson
    outname = fitted_tree_path + rawname + "_poisson"
    run_save_iqtree(alignment_file, "AA", "Poisson", outname1, true_tree_file, threads)
             

main()
