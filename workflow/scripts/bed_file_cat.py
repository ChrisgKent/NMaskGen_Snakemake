import csv


def main(bed_dict=snakemake.input, output=snakemake.output[0]):
    """
    Give this the bed file dictionary. Where the linerage is the key and the path top the file is the value

    This combines the files
    """
    mutates_pos_set = set()
    ref_name = set()

    for file_dir in bed_dict:
        with open(file_dir, mode="r") as file:
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
    with open(output, mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerows(bed_file)


if __name__ == "__main__":
    main()
