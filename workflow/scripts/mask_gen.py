from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
import csv


def main(
    bed_dir=snakemake.input["bed_dir"],
    ref_dir=snakemake.input["ref_dir"],
    output_dir=snakemake.output[0],
):
    """
    Given a list of positions and a referance sequence. This will mask all positions with an N.
    The mutated sequence is returned as a Seq object
    """
    # Read in the combined bed_file
    with open(bed_dir) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        positions = [x[1] for x in tsv_file]

    # Reads in the referance genomes and assigns the name and seq
    ref_it = SeqIO.parse(ref_dir, "fasta")
    ref_seqrecord = [x for x in ref_it]

    mut_seq = MutableSeq(ref_seqrecord[0].seq)

    for i in positions:
        mut_seq[int(i)] = "N"

    SeqIO.write(
        SeqRecord(mut_seq, id="NMask", description=""),
        output_dir,
        "fasta",
    )


if __name__ == "__main__":
    main()
