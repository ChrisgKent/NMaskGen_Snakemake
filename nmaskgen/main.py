import click
import pathlib
import multiprocessing

# Used in the insert loop
from itertools import compress

# Used for writing the .bed file
import csv

from click.core import F


from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo
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


def consen_gen(aligned_path, fasta_seq_name):
    """
    This takes the path to an aligned file and the name you want the consensus sequence to call called.
    Returns a SeqRecord containg the consensus sequence
    """
    align = AlignIO.read(aligned_path, "fasta")
    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.gap_consensus(threshold=0.5, ambiguous="N")

    consensus_seq_r = SeqRecord(Seq(consensus), id=fasta_seq_name, description="")
    return consensus_seq_r


def consen_gen2(aligned_path, fasta_seq_name, threshold=2 / 3):
    """ "
    This takes the path to an aligned file and the name you want the consensus sequence to call called.
    Returns a SeqRecord containg the consensus sequence

    This is differance from consen_gen. As this applies logic to which base is returned in a few ways

    Seq1:               ATGC--
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

    """
    align = AlignIO.read(aligned_path, "fasta")

    cat_seq = str()
    for seqrecord in align:
        cat_seq += str(seqrecord.seq)

    unique_bases = set(cat_seq)
    returned_base = [""] * align.get_alignment_length()
    for pos in range(align.get_alignment_length()):
        slice = ""
        for i in range(len(align)):
            slice += align[i][pos]

        consen_dict = {}
        for base in unique_bases:
            consen_dict[base] = slice.count(base)

        # Removes bases that confuse the analysis
        consen_dict.pop("N", None)
        consen_dict.pop("-", None)

        count_of_returned_base = max(consen_dict.values())
        max_base = {s for s in consen_dict if consen_dict[s] == count_of_returned_base}
        number_valid_base = sum(consen_dict.values())

        if number_valid_base == 0:
            returned_base[pos] = "N"
            # If there are no valid bases, then this slice must only contain "N" or "-".
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
    seq = SeqRecord(Seq("".join(returned_base)), id=fasta_seq_name, description="")
    return seq


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
    return dict(zip(x, group_list))


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
                fasta_seq_name=pango_lin + "_pseudo_consensus",
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

        ############################################################
        ## Starts to repair the left/right ends of the pre_repair ##
        ############################################################

        repair_test = AlignIO.read(pre_repair_consen_file, format="fasta")

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
            print(pango_lin, "has insert")

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

        SeqIO.write(
            SeqRecord(mutable_seq, id=pango_lin + "_pseudoref", description=""),
            output_pango_dir / pathlib.Path(pango_lin + "_pseudoref.fasta"),
            "fasta",
        )

        # Generates a .bed file, which contains all the changes between the p
        total_row_list = [""] * len(mutable_seq)
        diff_base_index = [""] * len(mutable_seq)

        for i in range(len(mutable_seq)):
            total_row_list[i] = [
                ref_name,
                i,
                i + 1,
                ref_sequence[i] + "to" + mutable_seq[i],
            ]
            diff_base_index[i] = ref_sequence[i] != mutable_seq[i]

        # Differant base?
        diff_bases = list(compress(total_row_list, diff_base_index))
        bed_dir = output_pango_dir / pathlib.Path(pango_lin + "_pseudoref.bed")
        with open(bed_dir, "w", newline="") as file:
            writer = csv.writer(file, delimiter="\t")
            writer.writerows(diff_bases)


if __name__ == "__main__":
    pass
