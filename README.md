#  R/Statest: statistical association test for continuous genotypes and phenotypes with R
 
Sunhee Kim and Chang-Yong Lee

R/Statest represents R scripts to perform Kendall’s and Pearson’s association test to assess the robustness under different genotype encodings by using simulated as well as real data sets.  

We proposed Kendall's test as a robust association test under different genotype encodings in a genome-wide association analysis of continuous traits. The R scripts provide the assessment of Kendall’s test and compares with that of Pearson's test in terms of the difference in p-values obtained by using different genotype encodings. We provide the R scripts together with real data sets in order for the readers to reproduce the results discussed in the manuscript.

**Loading the scripts**: copy gen_sim.R and stest.R from its GitHub repository
* gen_sim.R : generating simulated data
* stest.R : statistical hypothesis test (Kendall’s and Pearson’s test) and statistical analysis

**Install prerequisites**: install seqinr and Kendall packages 

     install.packages(“seqinr”)
     library(seqinr)
     install.packages(“Kendall”)
     library(Kendall)

     
## R-script tutorial
### Generating simulated data: gen_sim.R
To generate simulated data, run gen_sim.R with the following parameters.

     gen_sim(q, pi, snp, n, rho)
     #q: minor allele frequency
     #pi: the rate of the variation
     #snp : number of snp 
     #n: number of samples
     #rho: shuffling rate
     
(eg) gen_sim(0.25, 0.05, 20000, 300, 0.5)

**output** 
* sim_geno1.txt: simulated genotype data with encodings E_1={0,1,2} 
* sim_geno2.txt : simulated genotype data with encodings E_1={-1,0,1}
* sim_pheno.txt : simulated phenotype data   


### Performing statistical tests 
Association analysis. Simply load in genotype and phenotype data into stest.R as geno and pheno variables. Output will be the frequency distribution of ∆p_i and skewness and kurtosis of the distribution. We provide genotype and phenotype of real data, which are listed and described as follows

Once the data has been loaded into R, the stest.R for association test can be run with a single command:

Run the stest.R for association test with the following parameters

        stest(geno1, geno2, pheno, n_pca=4)
        #geno1, geno2: genotype data, pheno: phenotype data, n_pca: number of principal components

(eg1) stest(sim_geno1.txt, sim_geno2.txt, sim_pheno.txt, n_pca=4)    
(eg2) stest(real_geno1.txt, real_geno2.txt, real_pheno.txt, n_pca=4)  


**Output**
* xxx.jpg : the frequency distribution of ∆p_i    # xxx : trait name
* boxplots.jpg : box plots of skewness and kurtosis   

## Real data  
* real_geno1.txt: genotype data of 7,551 SNPs across 193 samples with encoding E_1={0,1,2} 
* real_geno2.txt: genotype data of 7,551 SNPs across 193 samples with encoding E_2={-1,0,1} 
* real_pheno.txt : phenotype data of 30 traits across 193 samples   
