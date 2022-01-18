To install
```
git clone https://github.com/ChrisgKent/NMaskGen_Snakemake
```
Activate / install the conda enviroment
```
cd NMaskGen_Snakemake
conda env create -f NMaskGen.yaml

conda activate nmaskgen

snakemake --cores all --config input_dir={input} 
```
