# DISCLAIMER
This is an in-house R script for genotype calling. This code was used to analyze the data described in Comparative assessment of SNP genotyping assays for challenging forensic samples utilizing ancient DNA methods (Staadig et. al).
# Genotype_Call
The Genotype_Call R script can be used to call genotypes from sequence data generated from CLC Genomics Workbench. An example of input file format is shown below. The R script was developed under R version 4.3.1.

The following threshold parameters are editable:
het_bal1_cut (Heterozygote balance)
het_bal2_cut (Heterozygote balance)
het_bal3_cut (Heterozygote balance)
hom_bal_cut (Homozygote balance)
Q_cut (Q-score)
cov_hom_cut (Coverage Homozygote)
cov_het_cut (Coverage Heterozygote)
cov_max (Maximum coverage)
F_R_cut (Forward/Reverse read count balance)

# Example
This is an example to run the script: Genotype_Call(het_bal1_cut=0.5,het_bal2_cut=0.7,het_bal3_cut=0.1,hom_bal_cut=0.95,Q_cut=25,cov_hom_cut=10,cov_het_cut=5,cov_max=5000,F_R_cut=0.19,indata="Input.txt",outfile="Output")

The input file must contain the following columns:
"Chromosome", "Region", "Allele", "Count", "Forward/reverse balance", "Average quality". Each marker should have four rows, one per nucelotide (A, C, G, T)

## License

This code is released under the MIT License (OSI-approved).  
See the `LICENSE` file for details.
