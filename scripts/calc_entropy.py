"""
Calculates both the AA-level entropy of a) all site preferences (ie sim parameters), and b) all simulated alignments (averaged across whole alignment) 
"""
import os 
import pyvolve
import numpy as np
g = pyvolve.Genetics()


def codon_freqs_to_aa_freqs(codonfreqs):
    '''
        Codon codon frequencies to amino-acid frequencies.
    '''
    cf = dict(zip(g.codons, codonfreqs))
    aa_freqs = []
    for family in g.genetic_code:
        total = 0.
        for syn in family:
            total += cf[syn]
        aa_freqs.append(total)
    aa_freqs = np.array(aa_freqs)
    assert((1. - np.sum(aa_freqs)) <= 1e-8), "\nImproperly converted codon to amino acid frequencies."
    return aa_freqs



def calculate_entropy(f):
    return -1. * np.sum ( f * np.log(f) )

all_names = ["LAC", "NP", "HA", "HIV"]

prefpath = "../simulations/preferences/"
simpath = "../simulations/alignments/"
outpath = "../results_analysis/csv_files/"
outfile_preferences = outpath + "preferences_entropy.csv"
outfile_simulations = outpath + "simulations_entropy.csv"


############# Entropy of PREFERENCES ##############
# print("Preferences entropy")
# outstring = "name,site,entropy\n"
# for name in all_names:
#     prefs = np.loadtxt(prefpath + name + "_prefs.csv", delimiter=",")
#     i = 0
#     for siteprefs in prefs:
#         sitefit = np.log(  siteprefs/np.sum(siteprefs)   ) ## renormalize DMS, my tolerance is a bit more stringent
# 
#         m = pyvolve.Model("mutsel", {"fitness": sitefit})
#         
#         codon_freqs = m.extract_state_freqs()
#         freqs    = codon_freqs_to_aa_freqs(codon_freqs)
#         h = str(calculate_entropy(freqs))
#         outstring += ",".join([name, str(i),h]) + "\n"
#         i+=1
# with open(outfile_preferences, "w") as f:
#     f.write(outstring.strip())




############# Entropy of ACTUAL SIMULATIONS ##############
print("Simulations entropy")
outstring = "name,tree,rep,total_entropy\n"
fastas = [x for x in os.listdir(simpath) if x.endswith("_AA.fasta")]  
for fasta in fastas:
    stuff = fasta.split("_") # NP_salichos_rep3_AA.fasta
    f = pyvolve.ReadFrequencies("amino_acid", file = simpath + fasta)
    h = calculate_entropy( f.compute_frequencies() )
    outline = stuff[0] + "," + stuff[1] + "," + stuff[2].replace("rep","") + "," + str(h)
    outstring += outline + "\n"
    
with open(outfile_simulations, "w") as f:
    f.write(outstring.strip())


    
for name in all_names:
    print(name)
    prefs = np.loadtxt(prefpath + name + "_prefs.csv", delimiter=",")
    i = 0
    for siteprefs in prefs:
        sitefit = np.log(  siteprefs/np.sum(siteprefs)   ) ## renormalize DMS, my tolerance is a bit more stringent

        m = pyvolve.Model("mutsel", {"fitness": sitefit})
        
        codon_freqs = m.extract_state_freqs()
        freqs    = codon_freqs_to_aa_freqs(codon_freqs)
        h = str(calculate_entropy(freqs))
        outstring += ",".join([name, str(i),h]) + "\n"
        i+=1
with open("../results_analysis/preferences_site_entropy.csv", "w") as f:
    f.write(outstring.strip())