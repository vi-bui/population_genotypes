####parse VCF file

#imports
#use pysam to parse VCF file
from itertools import count
from types import GeneratorType
from pysam import VariantFile
import numpy as np
from sklearn import decomposition
import pandas as pd

#asign the VCF file and panel with populations into a variable
vcf_file = "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
panel_file = "phase1_integrated_calls.20101123.ALL.panel"

genotype = []
samples = []
rsids = []
with VariantFile(vcf_file) as vcf_read:
    counter = 0
    for read in vcf_read:
        counter += 1
        if counter % 100 == 0: #sample snp across Chr22
            alleles = [read.samples[x].allele_indices for x in read.samples] #get the alleles
            samples = [sample for sample in read.samples] #get the sample names for a single SNP
            genotype.append(alleles)
            rsids.append(read.id)
        if counter % 4943 == 0:
            print(counter)
            print(f'{round(100 * counter / 494328)}%')

#add the population code to the samples
with open(panel_file) as file:
    annotations = {} #{sample:population}
    for point in file:
        point = point.strip().split('\t')
        annotations[point[0]] = point[1] #first column in panel is sample IDs and second is pop code

#print(genotype)
#print(len(genotype)) #get number of samples
print(rsids)
#get the alleles into a numpy array
#needed to perform PCA later
genotypes = np.array(genotype)
print(genotypes.shape) #100 is the first hundred SNPs. 2504 represents the number of samples in the samples. 2 is the number of chromosomes

matrix = np.count_nonzero(genotypes, axis=2) #only get the first 100 SNPs and 2504 samples into a matrix
print(matrix.shape)

#########################################
#run PCA on samples

#transpose the matrix 
matrix = matrix.T
print(matrix.shape)

#running PCA
pca = decomposition.PCA(n_components=2)
pca.fit(matrix)
print(pca.singular_values_)
pca_transformed = pca.transform(matrix)
print(pca_transformed.shape)

#producing dataframe
df = pd.DataFrame(matrix, columns=rsids, index=samples)
df["Population code"] = df.index.map(annotations)
df.to_csv("matrix.csv")