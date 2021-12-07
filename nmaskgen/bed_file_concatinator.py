import csv
import pathlib


def bed_concat(bed_dict, output):
    """
    Give this the bed file dictionary. Where the linerage is the key and the path top the file is the value

    This combines the files
    """
    mutates_pos_set = set()
    ref_name = set()

    for file_dir in bed_dict.keys():
        with open(bed_dict[file_dir], mode="r") as file:
            bedfile = csv.reader(file, delimiter="\t")

            for lines in bedfile:
                mutates_pos_set.add(int(lines[1]))
                ref_name.add(lines[0])

    mutates_pos_list = list(mutates_pos_set)
    mutates_pos_list.sort()

    bed_file = [""] * len(mutates_pos_list)
    for i in range(len(mutates_pos_list)):
        bed_file[i] = [
            list(ref_name)[0],
            mutates_pos_list[i],
            mutates_pos_list[i] + 1,
        ]

    # Writes a bed file to the output dir
    with open(pathlib.Path(output / "concat_bed.bed"), mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerows(bed_file)

    return mutates_pos_list
