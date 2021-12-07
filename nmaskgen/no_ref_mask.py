from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
import pathlib
from Bio.Seq import Seq


def no_ref_mask(pseudo_ref_dict, path):

    pseudo_list = [""] * len(pseudo_ref_dict)

    for i in range(len(list(pseudo_ref_dict))):
        seq = SeqIO.parse(pseudo_ref_dict[list(pseudo_ref_dict)[i]], "fasta")

        pseudo_list[i] = [x for x in seq][0]

    SeqIO.write(pseudo_list, path, format="fasta")

    pseudo_align = AlignIO.read(path, format="fasta")

    cat_seq = str()
    for seqrecord in pseudo_align:
        cat_seq += str(seqrecord.seq)

    # Finds all the bases contained within the alignments
    unique_bases = set(cat_seq)
    returned_base = [""] * pseudo_align.get_alignment_length()

    for pos in range(pseudo_align.get_alignment_length()):
        slice = ""
        for i in range(len(pseudo_align)):
            slice += pseudo_align[i][pos]

        consen_dict = {}

        # Counts now many times each base apears in each positional slice.
        # Saves results in a dict. With base as key and count as value
        for base in unique_bases:
            consen_dict[base] = slice.count(base)

        count_of_returned_base = max(consen_dict.values())
        max_base = {s for s in consen_dict if consen_dict[s] == count_of_returned_base}

        if count_of_returned_base == len(pseudo_align):
            returned_base[pos] = list(max_base)[0]
            # If there are no valid bases, then this slice must only contain "N" or "-".
        else:
            returned_base[pos] = "N"

    seq = SeqRecord(Seq("".join(returned_base)), id="nmask", description="")
    return seq
