# PCA on Genotypes

This repository consists of performing PCA and tSNE on the Phase 1 1000 Genomes Project data by samples and populations.

The VCF file and panel file containing the populations was downloaded from The International Genome Sample Resource using AWS by the following code

````
curl -O https://1000genomes.s3.amazonaws.com/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

curl -O https://1000genomes.s3.amazonaws.com/release/20110521/phase1_integrated_calls.20101123.ALL.panel


