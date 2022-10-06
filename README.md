To install
```
git clone https://github.com/ChrisgKent/NMaskGen_Snakemake
```
Activate / install the conda enviroment
```
cd NMaskGen_Snakemake
conda env create -f NMaskGen.yaml

conda activate nmaskgen
```
Running nmaskgen
```
snakemake --cores all --config input_dir={input} 
```
The example data
```
snakemake --cores all --config output_dir=example_data/new_results/ input_dir=example_data/sequences
```

Options

```input_dir```: Provides a directory, containing sub dirs for each varient

```ref_dir```: Provides a location of a fasta file containing the referance genome

```output_dir```: The directory the results will be saved to

