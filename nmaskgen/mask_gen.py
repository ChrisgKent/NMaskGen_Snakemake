from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord


def mask_gen(list_of_pos, ref):

    Mut_seq = MutableSeq(ref)

    for i in range(len(Mut_seq)):
        if i in list_of_pos:
            Mut_seq[i] = "N"

    return Mut_seq
