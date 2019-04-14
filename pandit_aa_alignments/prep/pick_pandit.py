from Bio import SeqIO
import random
import os

min_seq = 20
max_seq = 500
min_len = 100
max_len = 1000

choose  = 200
outdir = "use_pandit/"
indir = "pandit_fasta_aa/"
fastas = [x for x in os.listdir(indir) if x.endswith("fasta")]

contenders_sites = {}
contenders_seqs  = {}
contenders = []
for fasta in fastas:
    with open(indir) + fasta, "r") as f:
        recs = list(SeqIO.parse(f, "fasta"))
    nseq = len(recs)
    nsites = len(recs[0])
    
    seq_ok = nseq >= min_seq and nseq <= max_seq
    len_ok = nsites >= min_len and nsites <= max_len
    if seq_ok and len_ok:
        contenders.append(fasta)
        contenders_sites[fasta] = str(nsites)
        contenders_seqs[fasta] = str(nseq)
     

keepers = random.sample(list(range(len(contenders))), choose)

with open(outdir + "info.csv", "w") as f:
    f.write("name,nsites,nseq\n")
    for i in keepers:
        fasta = contenders[i]
        os.system("cp " + indir + fasta + " " + outdir)
        f.write(fasta.replace(".fasta", "") + "," + contenders_sites[fasta] + "," + contenders_seqs[fasta] + "\n")

    