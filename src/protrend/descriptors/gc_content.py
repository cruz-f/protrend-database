from Bio.SeqUtils import GC

from protrend.transform.functional_tfbs.tfbs import TFBSTransformer

def GCContent (motif):
    motif = TFBSTransformer.align_tfbs()
    descriptor = GC(motif)
    return descriptor