#  R/Statest: statistical association test for continuous genotypes and phenotypes with R
 
Sunhee Kim and Chang-Yong Lee

R/Statest represents R scripts to perform Kendall’s and Pearson’s association test to assess the robustness under different genotype encodings by using simulated as well as real data sets. 
We proposed Kendall's test as a robust association test under different genotype encodings in a genome-wide association analysis of continuous traits. The R scripts provide the assessment of Kendall’s test and compares with that of Pearson's test in terms of the difference in p-values obtained by using different genotype encodings. We provide the R scripts together with real data sets in order for the readers to reproduce the results discussed in the manuscript.

##Loading the scripts: copy gen_sim.R and stest.R from its GitHub repository
##Install prerequisites: install kendall packages
install.packages(“kendall”)
library(kendall)

##R-script tutorial
* Generating simulated data: gen_sim.R
To generate simulated data, run gen_sim.R with the following parameters.
q: minor allele frequency,   # pi: the rate of the variation,  # n: number of samples,	# rho: shuffling rate
gen_sim(q, pi, n, rho, sim_geno.txt, sim_pheno.txt)
output 
sim_geno.txt  # simulated genotype data
sim_pheno.txt  # simulated phenotype data

* Performing statistical tests 
Association analysis. Simply load in genotype and phenotype data into stest.R as geno and pheno variables. Output will be the frequency distribution of ∆p_i and skewness and kurtosis of the distribution. We provide genotype and phenotype of real data, which are listed and described as follows
Once the data has been loaded into R, the stest.R for association test can be run with a single command:
### Run the stest.R for association test with the following parameters
### geno: genotype data,	# pheno: phenotype data,	# n_pca: number of principal components
stest(geno1,geno2, pheno, n_pca=4)
Input 
geno: sim_geno1.txt, sim_geno2.txt (or real_geno1.txt, real_geno2.txt)
pheno: sim_pheno.txt (or, real_pheno.txt) 

Output
tname.jpg : the frequency distribution of ∆p_i   # tname : trait name 
boxplot.jpg : box plots of skewness and kurtosis

Real data
real_geno1.txt, real_geno2.txt : genotype data of 7,551 SNPs across 193 samples with encoding E_1={0,1,2} and E_2={-1,0,1} 
real_pheno.txt : phenotype data of 30 traits across 193 samples


