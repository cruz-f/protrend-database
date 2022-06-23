from Bio.motifs.matrix import PositionSpecificScoringMatrix
from Bio.Alphabet import IUPAC

from protrend.transform.functional_tfbs.tfbs import TFBSTransformer

def PSSM_BP (motif):
    motif = TFBSTransformer.align_tfbs()
    descriptor = PositionSpecificScoringMatrix(alphabet=IUPAC, values=motif)
    return descriptor


# ou

from Bio import motifs

def PSSM (motif):
    motif = TFBSTransformer.align_tfbs()
    m = motifs.create(motif)
    pwm = m.counts.normalize(pseudocounts=0.5)
    descriptor2 = pwm.log_odds()
    return descriptor2