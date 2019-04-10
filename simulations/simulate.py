"""
    SJS
    Simulate alignments according to DMS parameters
"""
import os
import sys
import numpy as np
import pyvolve
from copy import deepcopy 
# sys.path.append("../")
# from sw_functions import *


translation_table = {
                    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
                    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                    'TAC':'Y', 'TAT':'Y', 'TGC':'C', 'TGT':'C', 
                    'TGG':'W'
                    }
 


names = [ "HIV", "Gal4", "LAC", "HA", "NP"]
mudict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}

simrep = sys.argv[1] ## 1-20

treepath = "true_trees/"
prefpath = "preferences/"
simpath  = "alignments/"
trees = [x for x in os.listdir(treepath) if x.endswith("bl0.03_rep1.tree")]


all_partitions = {}

for name in names: ### 5
    prefs = np.loadtxt(prefpath + name + "_prefs.csv", delimiter=",")
    partitions = []
    for siteprefs in prefs:
        sitefit = np.log(  siteprefs/np.sum(siteprefs)   ) ## renormalize, my tolerance is a bit more stringent
        m = pyvolve.Model("mutsel", {"fitness": sitefit, "mu": mudict})
        p = pyvolve.Partition(models = m, size = 1)
        partitions.append(p)
    all_partitions[name] = partitions

for tree in trees:  ### 3 for now
    with open(treepath + tree, "r") as f:
        treestring = f.read().strip()

    pytree = pyvolve.read_tree(tree = treestring)
    treename = tree.split(".tree")[0]
    
    for name in all_partitions:
        
        outname = simpath + name + "_" + treename + "_rep" + str(simrep) + "_CODON.fasta" ##### CODON
        outname1 = outname.replace("_CODON", "_DNA")               ##### NUCLEOTIDE (same as codon but naming is convenient downstream)
        outname2 = outname.replace("_CODON", "_AA")                ###### AMINO ACID
        outname3 = outname2.replace(".fasta", ".dat")              ###### AMINO ACID

        if os.path.exists(outname):
            continue

        e = pyvolve.Evolver(partitions = deepcopy(all_partitions[name]), tree = deepcopy(pytree))
        e(seqfile = outname, seqfmt = "fasta", ratefile = None, infofile = None)

        nucseqs = e.get_sequences()
        aaseqs = {}
        for id in nucseqs:
            seq = nucseqs[id]
            protein = ""    
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3]
                protein += translation_table[codon]
            aaseqs[id] = protein

        os.system("cp " + outname + " " + outname1)

        aa_out = ""
        for record in aaseqs:
            aa_out += ">" + str(record) + "\n" + str(aaseqs[record]) + "\n"

        with open(outname2, "w") as f:
            f.write(aa_out)

        with open(outname3, "w") as f:
            f.write(aa_out)
            f.write("\n" + treestring) 
