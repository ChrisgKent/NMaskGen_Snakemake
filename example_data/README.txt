Within ~/sequences there are 4 example variants, which are all based on the reference sequence.  

Sequences within each variant all share the base mutations (random SNPS) as well as each sample having a random number of point mutations. Additionally, a random number of bases (between 0-30) are removed from each end to mimic the poor coverage at the 3' and 5' end.

Input: 
The input sequences need to have a similar directory structure. 

This example was generated with:
```
snakemake --cores all --config output_dir=example_data/results/ input_dir=example_data/sequences
```


