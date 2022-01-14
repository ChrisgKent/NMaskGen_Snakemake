from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main(
    aligned_path=snakemake.input[0],
    output=snakemake.output[0],
    seq_name=snakemake.params["clade"],
    threshold=snakemake.params["threshold"],
    ignore_n=True,
):
    """
    This takes the path to an aligned file and the name you want the consensus sequence to call called.
    Returns a SeqRecord containg the consensus sequence

    This is differance from consen_gen. As this applies logic to which base is returned in a few ways

    Seq1:               ATGCA-
    Seq2:               ATCC--
    Seq3:               N-----
    Seq4:               N-----

    classic_consen:     N-----
    this_consen:        ATNCNN


    A simple consen_gen would have returned "N---" as "-" is the most common base in each position.
        Expect for the first position where the equal number would trigger an ambiquity char

    This consen discards the gaps (or "N"), enabling genomes with poly N track to be used.
        It can also detect vairability within the valid bases
        If no valid bases are detected it will return "N"

    However it does require, two sequences to generate a base. To prevent point mutations being represented

    """
    align = AlignIO.read(aligned_path, "fasta")

    cat_seq = str()
    for seqrecord in align:
        cat_seq += str(seqrecord.seq)

    # Finds all the bases contained within the alignments
    unique_bases = set(cat_seq)
    returned_base = [""] * align.get_alignment_length()
    # For each base position, it takes a slice
    for pos in range(align.get_alignment_length()):
        slice = ""
        for i in range(len(align)):
            slice += align[i][pos]

        consen_dict = {}
        # Counts now many times each base apears in each positional slice.
        # Saves results in a dict. With base as key and count as value
        for base in unique_bases:
            consen_dict[base] = slice.count(base)

        # Removes bases that confuse the analysis. However these are consideered later
        if ignore_n:
            consen_dict.pop("N", None)

        count_of_returned_base = max(consen_dict.values())
        max_base = {s for s in consen_dict if consen_dict[s] == count_of_returned_base}
        number_max_base = consen_dict[list(max_base)[0]]

        number_valid_base = sum(consen_dict.values())

        if number_valid_base == 0:
            returned_base[pos] = "N"
            # If there are no valid bases, then this slice must only contain "N" or "-".
        elif len(max_base) == 1 and number_max_base == 1:
            returned_base[pos] = "N"
        elif len(max_base) == 1:
            returned_base[pos] = list(max_base)[0]
            # If there is only one key, which has the highest value. Then it is assigned to the position
        elif len(max_base) != 1:
            returned_base[pos] = "N"
            # If two valid bases both have the same ocourance. Then an "N" is asigned

            # If the assigned base doesn't meet the threshold. Then it is replaced with an "N"
        if (
            0 != number_valid_base
            and count_of_returned_base / number_valid_base < threshold
        ):
            returned_base[pos] = "N"
    # The list is joined into a string
    seq = SeqRecord(Seq("".join(returned_base)), id=seq_name, description="")

    # Write the output file
    SeqIO.write(
        seq,
        output,
        "fasta",
    )


if __name__ == "__main__":
    main()
