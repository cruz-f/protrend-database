from Bio.motifs.matrix import PositionWeightMatrix
from Bio.Alphabet import IUPAC

from protrend.transform.functional_tfbs.tfbs import TFBSTransformer

def PWM_BP (motif):
    motif = TFBSTransformer.align_tfbs()
    descriptor = PositionWeightMatrix(alphabet=IUPAC, values=motif)
    return descriptor


# ou

from Bio import motifs

def PWM (motif):
    motif = TFBSTransformer.align_tfbs()
    m = motifs.create(motif)
    descriptor1 = m.counts.normalize(pseudocounts=0.5)
return descriptor1