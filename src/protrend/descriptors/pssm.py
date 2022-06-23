from Bio.motifs.matrix import PositionSpecificScoringMatrix
from Bio.Alphabet import IUPAC

from protrend.transform.functional_tfbs.tfbs import TFBSTransformer

def PSSM (motif):
    motif = TFBSTransformer.align_tfbs()
    descriptor = PositionSpecificScoringMatrix(alphabet=IUPAC, values=motif)
    return descriptor