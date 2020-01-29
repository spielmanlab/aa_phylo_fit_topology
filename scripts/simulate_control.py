"""
    SJS
    Simulate alignments using a _standard model_ of protein evolution: WAG+I+G where p(I) = 0.05 and \Gamma ASRV with 10 categories, \alpha = 0.8
"""
import os
import sys
import numpy as np
import pyvolve
from copy import deepcopy 

treepath = "../simulations/true_trees/" ## since simulating nucleotide level but want BL to describe, more or less, protein divergence
prefpath = "../simulations/preferences/"
simpath  = "../simulations/alignments/wag_control/"
trees = [x for x in os.listdir(treepath) if x.endswith("_resolved.tree")] ## 8 trees

for name in ["LAC", "NP", "HA", "HIV"]:
    print("Setting up model for for", name)
    prefs = np.loadtxt(prefpath + name + "_prefs.csv", delimiter=",") ## just to get length
    
    m = pyvolve.Model("WAG", alpha = 0.8, num_categories = 10, pinv = 0.05)
    partitions = pyvolve.Partition(models = m, size = 1)#len(prefs))
    
    for tree in trees:
        with open(treepath + tree, "r") as f:
            treestring = f.read().strip()
        pytree = pyvolve.read_tree(tree = treestring)
        treename = tree.replace("_resolved.tree","")
        print("Simulating along", treename)
    
        for i in range(1, 21):
            outname = simpath + name + "_" + treename + "_rep" + str(i) + "_AA.fasta"
    
            if os.path.exists(outname):
                continue

            e = pyvolve.Evolver(partitions = deepcopy(partitions), tree = deepcopy(pytree))
            e(seqfile = outname, seqfmt = "fasta", ratefile = None, infofile = None)

