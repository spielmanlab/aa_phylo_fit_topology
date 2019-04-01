import os
import sys
import pyvolve
from sw_functions import *
import sys
import numpy as np
 
"""
    This script calculates an ANALOGY for the protein exchangability model.
    
    Exchangeability model can be decomposed as Qij = Rij*F_j, where Rij is an exchangeability ("rate") and F_j is the target frequency.
    For a given alignment, swMutSel has inferred site-wise fitness and global mutation rates. We can use these to get site-wise Rij. 
    For making the models we will put 0.05's all around and just do +F for inference.   
    
    Save Q matrix (with 0 diagonal) in text format.
    Save PAML model (decomposed rates + dummy 0.05 freqs from the HB model) to text.
"""


ZERO    = 1e-8
G            = pyvolve.Genetics()
GENETIC_CODE = G.genetic_code
CODONS       = G.codons
SWPATH     = "simulations/swmutsel_output/"
OUTPATH    = "hb_models/"
PAML_AA    = "ARNDCQEGHILKMFPSTWYV"
PYVOLVE_AA = "ACDEFGHIKLMNPQRSTVWY"


def codon_freqs_to_aa_freqs(codonfreqs):
    '''
        Codon codon frequencies to amino-acid frequencies.
    '''
    g = pyvolve.Genetics()
    cf = dict(zip(g.codons, codonfreqs))
    aa_freqs = []
    for family in g.genetic_code:
        total = 0.
        for syn in family:
            total += cf[syn]
        aa_freqs.append(total)
    aa_freqs = np.array(aa_freqs)
    assert((1. - np.sum(aa_freqs)) <= ZERO), "\nImproperly converted codon to amino acid frequencies."
    return aa_freqs



def main():
    swfiles = [x for x in os.listdir(SWPATH) if x.endswith("MLE.txt")]
    
    paml_frequencies = "\n" +  " ".join(["0.05"]*20)
    for rawname in swfiles:
        name = rawname.split("_")[0]
        outfilepaml = OUTPATH + name + "_HB.paml" 
        print(name)
    
        with open(SWPATH + rawname, "r") as f:
            swstuff = f.readlines()
        mu = parse_swmutsel_mutation(swstuff)
        all_fitness = parse_swmutsel_fitness(swstuff)
        
        ## Average codon rates
        mean_codon_matrix = np.zeros([61, 61])
        for site_fitness in all_fitness:
            params = {"fitness": site_fitness, "mu": mu}
            m = pyvolve.Model("mutsel", params)
            mean_codon_matrix += m.extract_rate_matrix() ## 61x61
        mean_codon_matrix /= len(all_fitness)
        
        
        ## Convert to amino acid model assuming equal frequencies
        aa_matrix = np.zeros([20,20])
        for i in range(20):
            source_codons = GENETIC_CODE[i]
            for j in range(i+1, 20):
                target_codons = GENETIC_CODE[j]
                
                aa_rate = 0.
                x = 0
                for s in source_codons:
                    source_index = CODONS.index(s)
                    for t in target_codons:
                        target_index = CODONS.index(t)
                        this_rate = mean_codon_matrix[source_index][target_index]
                        assert(abs(this_rate) == this_rate), "rate is negative somehow FREAK OUT" 
                        aa_rate += this_rate
                        x += 1
                aa_rate /= x
                
                aa_matrix[i][j] = aa_rate
                aa_matrix[j][i] = aa_rate
        for i in range(20):
            aa_matrix[i][i] = -1. * np.sum(aa_matrix[i])  
            
        ## Decompose to get rate matrix (not q matrix)
        decomped_aa_matrix = np.dot(aa_matrix, np.diag(1.0 / np.repeat(0.05, 20)))  ## frequencies will be equal based on how calculations above were done
        assert(np.allclose(decomped_aa_matrix, decomped_aa_matrix.T)), "\n decomp'd rate matrix not symmetric"

        ## Save to PAML format
        paml = ""
        full = np.zeros([20,20])
        for s in range(20):
            for t in range(20):
                if s==t:
                    continue
                saa = PAML_AA[s] ## paml
                taa = PAML_AA[t] ## paml
                full[s][t] = decomped_aa_matrix[ PYVOLVE_AA.index(saa) ][ PYVOLVE_AA.index(taa) ]
        assert(np.allclose(full, full.T)), "\n paml converting rate matrix not symmetric"
        for i in range(20):
            x = 0
            row = ""
            while x < i:
                row += str(full[x][i]) + " " 
                x+=1
            paml += row.strip() + "\n"

        paml += paml_frequencies
        with open(outfilepaml, "w") as f:
            f.write(paml.strip())

main()

   
