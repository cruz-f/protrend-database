from Bio import motifs
from Bio.Seq import Seq

def MotifLogo (aligned_seqs):
    seqs = []
    for seq in aligned_seqs:
        seqs.append(Seq(seq))
    m = motifs.create(seqs, alphabet='ATGC-')
    descriptor = m.weblogo("mymotif.png")
    return descriptor