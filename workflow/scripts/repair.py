from Bio import AlignIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from itertools import compress


def continous_func(x):
    """
    This take a list of intergers, and group them into continous runs.
    Returns a dictionry that contains the interger (key) and the group (value)
    """
    group = 0
    group_list = [""] * len(x)
    for insert_positions in range(len(x)):
        if insert_positions == 0:
            group_list[insert_positions] = group
        elif x[insert_positions] - x[insert_positions - 1] == 1:
            group_list[insert_positions] = group
        else:
            group += 1
            group_list[insert_positions] = group
    return dict(zip(x, group_list))


def main(
    aligned_psuedo_path=snakemake.input[0],
    output=snakemake.output[0],
    seq_name=snakemake.params["seq_name"],
):
    repair_test = AlignIO.read(aligned_psuedo_path, format="fasta")
    mutable_seq = MutableSeq(repair_test[0].seq)

    # Starting from 0 this iterates through, and if the sequences are differant it replaces with bases from the ref
    ## repair-test[1] is the ref seq

    for pos in range(repair_test.get_alignment_length()):
        agreement = set(repair_test[:, pos])
        if len(agreement) != 2:
            break
        mutable_seq[pos] = repair_test[1][pos]

    # Starts from the right.
    for pos in range(repair_test.get_alignment_length())[::-1]:
        agreement = set(repair_test[:, pos])
        if len(agreement) != 2:
            break
        mutable_seq[pos] = repair_test[1][pos]

    # Detects delections within the psudeogenome
    del_indexes = [i for i, j in enumerate(mutable_seq) if j == "-"]
    for i in del_indexes:
        mutable_seq[i] = "N"

    # Detects insertions into the pseudogenome
    # If there has been an insertion into the ref genome
    if repair_test[1].seq.count("-") > 0:
        # This contains the indexes of all locations of "-" in the aligned ref genome
        indexes = [i for i, j in enumerate(repair_test[1].seq) if j == "-"]
        # This function return if the insertions are adjasent. "indexes" key is the index, the value is the group
        indexes = continous_func(indexes)

        # Solves each group
        for insert_group in set(indexes.values()):
            # Finds all indexes that belong to the group "insert_group"
            for_index = [x == insert_group for x in list(indexes.values())]
            values = list(compress(list(indexes.keys()), for_index))

            # Replaces all of the bases the deleted with an *
            for i in values:
                mutable_seq[i] = "*"

            # Replaces the base -1 with an "N"
            mutable_seq[min(values) - 1] = "N"

        # Once all the bases have been marked, they can then be deleted without messing up the index system
        for i in range(mutable_seq.count("*")):
            mutable_seq.remove("*")

    seq = SeqRecord(mutable_seq, id=seq_name, description="")
    # Write the output file
    SeqIO.write(
        seq,
        output,
        "fasta",
    )


if __name__ == "__main__":
    main()
