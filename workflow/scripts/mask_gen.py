from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord


def mask_gen(list_of_pos, ref):
    """
    Given a list of positions and a referance sequence. This will mask all positions with an N.
    The mutated sequence is returned as a Seq object
    """
    Mut_seq = MutableSeq(ref)

    for i in range(len(Mut_seq)):
        if i in list_of_pos:
            Mut_seq[i] = "N"

    return Mut_seq
