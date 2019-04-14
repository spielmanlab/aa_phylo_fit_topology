from Bio import SeqIO
import os
from urllib.request import urlopen

tree_url = "https://www.ebi.ac.uk/goldman-srv/pandit/pandit.cgi?action=browse&fam=REPLACE&field=rtp-rph"

def clean_string(st):
    return st.replace("-","").replace("/","").replace("_","").strip()


def obtain_source(search_url, type):

    response = urlopen(search_url)
    if type == "lines":
        page_source = [x.decode(response.headers.get_content_charset()) for x in response.readlines()]
    elif type == "full":
        page_source = response.read().decode(response.headers.get_content_charset())
        
    return page_source



fastas = [x for x in os.listdir(".") if x.endswith("fasta")]

for fasta in fastas:
    print(fasta)
    
    records = list(SeqIO.parse(fasta, "fasta"))
    with open(fasta, "w") as f:
        for record in records:
            f.write(">" + clean_string(record.id) + "\n" + str(record.seq) + "\n")
           
    name = fasta.replace(".fasta", "")

    url = tree_url.replace("REPLACE",name)
    page_source = obtain_source(url, "lines")

    for line in page_source:
        if line.startswith('<td align="left" nowrap>'):
            ts = line.replace('<td align="left" nowrap>', "").replace("</td>",'')
    
    treefile = name + ".tree"
    datfile  = name + ".dat"
    with open(treefile, "w") as f:    
        f.write(clean_string(ts))
    
    os.system("cat "+ fasta + " " + treefile + " > " + datfile)