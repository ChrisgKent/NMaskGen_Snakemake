import click
import pathlib
import multiprocessing

# Used for writing the .bed file
import csv

from click.core import F

# Imports consen_gen2()
from nmaskgen.consensus_gen import *
from nmaskgen.mask_gen import mask_gen

# Imports repair() and continous_func()
from nmaskgen.repair import *

# Imports bed_concat
from nmaskgen.bed_file_concatinator import *


from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord


def check_or_create_outpath(path, force: bool = False):
    """
    Check for an existing output dir, require --force to overwrite.
    Create dir if required, return Path obj.
    """
    if path.exists() and not force:
        raise IOError("Directory exists add --force to overwrite")

    path.mkdir(exist_ok=True)
    return path


@click.command()
@click.argument(
    "input",
    nargs=1,
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path
    ),
)
@click.argument(
    "ref_genome",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        path_type=pathlib.Path,
    ),
)
@click.option(
    "--output",
    type=click.Path(
        file_okay=False,
        dir_okay=True,
        writable=True,
        path_type=pathlib.Path,
    ),
    default="output",
)
@click.option(
    "--cores",
    type=click.IntRange(min=1, max=multiprocessing.cpu_count()),
    default=multiprocessing.cpu_count(),
)
def main(input, ref_genome, output=None, cores=None):
    check_or_create_outpath(output, force=True)

    # Creates a dict containing .bed file locations
    bed_file_dict = {}
    pseudo_ref_dict = {}

    # Reads in the referance sequence
    ref_fasta = SeqIO.parse(ref_genome, "fasta")

    pango_dirs = {}
    for child in input.iterdir():
        files = []
        if child.is_dir():
            for file in child.iterdir():
                if file.suffix == ".fasta":
                    files.append(list(SeqIO.parse(file, "fasta"))[0])
            pango_dirs[child.parts[-1]] = files

    for pango_lin in pango_dirs:
        # Sets the path to the dir to use as output in this interation of the loop
        output_pango_dir = pathlib.Path(output / pango_lin)
        # Creates the dir if it doesn't exist
        if not output_pango_dir.is_dir():
            output_pango_dir.mkdir()

        # Creates a tmp subdirectory, to dump intermediate files.
        output_pango_dir_tmp = output_pango_dir / pathlib.Path("tmp")
        if not output_pango_dir_tmp.is_dir():
            output_pango_dir_tmp.mkdir()

        SeqIO.write(
            pango_dirs[pango_lin],
            output_pango_dir / pathlib.Path(pango_lin + "_genomes.fasta"),
            "fasta",
        )

        # To save my sanity, MSA is only run if the align file doesn't already exist
        msa_file = pango_lin + "_pseudo_msa.fasta"
        msa_path = pathlib.Path(output_pango_dir_tmp / msa_file)
        if msa_path.is_file():
            print(msa_path, "exists")
        else:
            clustalomega_cline = ClustalOmegaCommandline(
                infile=output_pango_dir / pathlib.Path(pango_lin + "_genomes.fasta"),
                outfile=output_pango_dir_tmp
                / pathlib.Path(pango_lin + "_pseudo_msa.fasta"),
                verbose=True,
                percentid=True,
                threads=cores,
            )
            msa = clustalomega_cline()

        # This creates an aligned object, and then a consensus sequence
        pseudo_consensus = consen_gen2(
            output_pango_dir_tmp
            / pathlib.Path(
                pango_lin + "_pseudo_msa.fasta",
            ),
            fasta_seq_name=pango_lin + "_pseudo_consensus",
        )

        # Combines the consensus and the ref. Writes to tmp
        for seq_record in ref_fasta:
            ref_sequence = seq_record.seq
            ref_name = seq_record.id
        seq_ref = [pseudo_consensus, SeqRecord(ref_sequence, id=ref_name)]
        seq_ref_file = pango_lin + "_seqref.fasta"
        SeqIO.write(
            seq_ref,
            pathlib.Path(output_pango_dir_tmp / seq_ref_file),
            "fasta",
        )

        pre_repair_consen_file = output_pango_dir_tmp / pathlib.Path(
            pango_lin + "_repair_msa.fasta"
        )

        if pre_repair_consen_file.is_file():
            print(
                pre_repair_consen_file,
                "exists",
            )
        else:
            # Aligns the pseudo genome and the ref sequence
            clustalomega_cline = ClustalOmegaCommandline(
                infile=pathlib.Path(output_pango_dir_tmp / seq_ref_file),
                outfile=pre_repair_consen_file,
                verbose=True,
                percentid=True,
                threads=cores,
            )
            msa = clustalomega_cline()

        repaired_genome = repair(pre_repair_consen_file, pango_lin=pango_lin)

        SeqIO.write(
            SeqRecord(repaired_genome, id=pango_lin + "_pseudoref", description=""),
            output_pango_dir / pathlib.Path(pango_lin + "_pseudoref.fasta"),
            "fasta",
        )

        pseudo_ref_dict[pango_lin] = output_pango_dir / pathlib.Path(
            pango_lin + "_pseudoref.fasta"
        )

        # Generates a .bed file, which contains all the changes between the p
        total_row_list = [""] * len(repaired_genome)
        diff_base_index = [""] * len(repaired_genome)

        for i in range(len(repaired_genome)):
            total_row_list[i] = [
                ref_name,
                i,
                i + 1,
                ref_sequence[i] + "to" + repaired_genome[i],
            ]
            diff_base_index[i] = ref_sequence[i] != repaired_genome[i]

        # Differant base?
        diff_bases = list(compress(total_row_list, diff_base_index))
        bed_dir = output_pango_dir / pathlib.Path(pango_lin + "_pseudoref.bed")

        with open(bed_dir, mode="w", newline="") as file:
            writer = csv.writer(file, delimiter="\t")
            writer.writerows(diff_bases)

        bed_file_dict[pango_lin] = bed_dir

    concat_bed = bed_concat(bed_file_dict, output=output)

    pseudo_ref_dict[ref_name] = ref_genome

    changes = [""] * (len(concat_bed))

    psudo_test = [""] * len(pseudo_ref_dict)
    i = 0
    for lin in pseudo_ref_dict.keys():
        x = SeqIO.parse(pseudo_ref_dict[lin], format="fasta")
        for seqr in x:
            psudo_test[i] = seqr

        i += 1

    for i in range(len(concat_bed)):
        combined_bases = ""
        for seq in psudo_test:
            combined_bases += seq.seq[concat_bed[i]]

        changes[i] = combined_bases

    # Defining the headings
    y = [x.id for x in psudo_test]
    y = ["Position"] + y

    # Defining the rows
    x = [list(x) for x in changes]
    for i in range(len(concat_bed)):
        x[i] = [concat_bed[i]] + x[i]

    with open(output / pathlib.Path("base_changes.tsv"), mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerows([y] + x)

    # Generate an NMasked_fasta
    mask = mask_gen(concat_bed, ref_sequence)

    SeqIO.write(
        SeqRecord(mask, id="Nmask", description=""),
        output / pathlib.Path("NMask.fasta"),
        "fasta",
    )


if __name__ == "__main__":

    main()
