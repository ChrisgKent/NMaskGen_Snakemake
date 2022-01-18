import pathlib
from Bio import SeqIO

def main(input = snakemake.input[0],output = snakemake.output[0]):
    
    files = []
    for file in pathlib.Path(input).iterdir():
        if file.suffix == ".fasta":
            files.append(list(SeqIO.parse(file, "fasta"))[0])
    
    SeqIO.write(
        files,
        output,
        "fasta",
        )

if __name__ == "__main__":
    main()