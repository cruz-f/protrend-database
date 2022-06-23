from Bio.motifs.matrix import PositionWeightMatrix
from Bio.Alphabet import IUPAC

from protrend.transform.functional_tfbs.tfbs import TFBSTransformer

def PSSM (motif):
    motif = TFBSTransformer.align_tfbs()
    descriptor = PositionWeightMatrix(alphabet=IUPAC, values=motif)
    return descriptor