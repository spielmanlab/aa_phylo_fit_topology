import os
fastas = [x for x in os.listdir(".") if x.endswith("fasta")]

outfile = "pandit_info.csv"
outstring = "name,nsites,ntaxa,percent_ambig\n"
for fasta in fastas:
    nambig = 0
    ntaxa = 0
    nsites = 0
    with open(fasta, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith(">"):
            ntaxa += 1
        else:
            nambig += line.count("-") + line.count("X") + line.count("B") + line.count("J") + line.count("Z")
            if nsites == 0:
                nsites = len(line.strip())
    
    final = nambig/(nsites*ntaxa)
    
    assert (nambig != 0)
    assert (nsites != 0)
    assert (ntaxa != 0)
    
    outstring += ",".join( [fasta.replace(".fasta", ""), str(nsites), str(ntaxa), str(final) ]) + "\n"

with open(outfile,"w") as f:
    f.write(outstring.strip())

