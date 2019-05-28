"""
    SJS
    Simulate alignments according to DMS parameters
"""
import os
import sys
import numpy as np
import pyvolve
from copy import deepcopy 


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
 


mudict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}

name   = sys.argv[1]
simrep = sys.argv[2]

treepath = "true_trees/" ## since simulating nucleotide level but want BL to describe, more or less, protein divergence
prefpath = "preferences/"
simpath  = "alignments/"
trees = [x for x in os.listdir(treepath) if x.endswith("_resolved.tree")] ## 8 trees
print(trees)
prefs = np.loadtxt(prefpath + name + "_prefs.csv", delimiter=",")
partitions = []

print("Loading parameters")
for siteprefs in prefs:
    sitefit = np.log(  siteprefs/np.sum(siteprefs)   ) ## renormalize, my tolerance is a bit more stringent
    m = pyvolve.Model("mutsel", {"fitness": sitefit, "mu": mudict})
    p = pyvolve.Partition(models = m, size = 1)
    partitions.append(p)


for tree in trees:
    with open(treepath + tree, "r") as f:
        treestring = f.read().strip()
    pytree = pyvolve.read_tree(tree = treestring)#, scale_tree=3)
    treename = tree.replace("_resolved.tree","")
    print("Simulating along", treename)
    
    #outname = simpath + name + "_" + treename + "_rep" + str(simrep) + "_CODON.fasta" ##### CODON
    #outname1 = outname.replace("_CODON", "_DNA")               ##### NUCLEOTIDE (same as codon but naming is convenient downstream)
    outname = simpath + name + "_" + treename + "_rep" + str(simrep) + "_AA.fasta"
    outname_dat = outname.replace(".fasta", ".dat")    
    
    if os.path.exists(outname):
        continue
    #else:
    #    print(outname)
    #continue
    e = pyvolve.Evolver(partitions = deepcopy(partitions), tree = deepcopy(pytree))
    e(seqfile = outname.replace("AA", "CODON"), seqfmt = "fasta", ratefile = None, infofile = None)

    nucseqs = e.get_sequences()
    aaseqs = {}
    for id in nucseqs:
        seq = nucseqs[id]
        protein = ""    
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += translation_table[codon]
        aaseqs[id] = protein

    aa_out = ""
    for record in aaseqs:
        aa_out += ">" + str(record) + "\n" + str(aaseqs[record]) + "\n"

    with open(outname, "w") as f:
        f.write(aa_out)

    with open(outname_dat, "w") as f:
        f.write(aa_out)
        f.write("\n" + treestring) 
