from Bio import SeqIO
from itertools import compress
import csv


def main(
    input_seq=snakemake.input["seq"],
    ref_seq=snakemake.input["ref"],
    bed_output=snakemake.output[0],
):

    # Reads in the repaired genome, and assigns the seq
    repaired_it = SeqIO.parse(input_seq, "fasta")
    repaired_seqrecord = [x for x in repaired_it]
    repaired_genome = repaired_seqrecord[0].seq

    # Reads in the referance genomes and assigns the name and seq
    ref_it = SeqIO.parse(ref_seq, "fasta")
    ref_seqrecord = [x for x in ref_it]
    ref_name = ref_seqrecord[0].id
    ref_genome = ref_seqrecord[0].seq

    # Generates a .bed file, which contains all the changes between the p
    total_row_list = [""] * len(repaired_genome)
    diff_base_index = [""] * len(repaired_genome)

    # Iterates through all bases.
    # diff_base_index contains a boolen where True means differance in bases
    for i in range(len(repaired_genome)):
        total_row_list[i] = [
            ref_name,
            i,
            i + 1,
            ref_genome[i] + "to" + repaired_genome[i],
        ]
        diff_base_index[i] = ref_genome[i] != repaired_genome[i]

    # Selects the positions in which there is differance.
    diff_bases = list(compress(total_row_list, diff_base_index))

    with open(bed_output, mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerows(diff_bases)


if __name__ == "__main__":
    main()
