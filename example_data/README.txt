Within ~/sequences there are 4 example variants, which are all based on the reference sequence.  

Sequences within each variant all share the base mutations (random SNPS) as well as each sample having a random number of point mutations. Additionally, a random number of bases (between 0-30) are removed from each end to mimic the poor coverage at the 3' and 5' end.

Input: 
The input sequences need to have a similar directory structure. 

This example was generated with:
poetry run nmaskgen --output example_data/output --consen_thresh 0.8 example_data/sequences example_data/example_ref.fasta

--consen_thresh is the value in which the base will be masked with an N. 
	For example; if 6/7 bases are conserved. The threshold is 0.86, which is greater than 0.8.
	Hence the position isn't masked.  
	If two bases are different (either the same or different SNPs), then the threshold of 5/7 (0.71) is less than the threshold and the position is masked.

Output:
Each variants has its own output containing;
	x_genomes.fasta: Contains all the input genomes
	x_pseudoref.bed:	Contains all positions of variation relative to the ref genome
	x_pseudoref.fasta: Contains the nmasked sequence for this pseudolin
	tmp/
	x_pseudo_msa:	Contains all MSA of all the raw genomes
	x_repair_msa:	Contains the MSA of the consensus genome and the reference
	x_seqref:	The pre-aligment file with the consensus genome and the reference

base_changes.tsv: 	Shows all the mutations, and the source
Concat_bed.bed:		The bed file for the nmask
nmask:			The combined mask	


