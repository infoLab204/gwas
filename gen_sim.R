if(!require(seqinr)) install.packages("seqinr")
library(seqinr)

gen_sim <- function(q, pi,snp,n, rho)
{
    
    geno_type1 <- matrix(data=NA, nrow=snp, ncol=n)
    start_geno <- matrix(NA, nrow=1, ncol=n)
    start_geno[1,] <- sample(0:2,n, replace=TRUE,prob=c((1-q)*(1-q),2*q*(1-q),q*q))

    pheno <- matrix(data=NA, nrow=n, ncol=1)
    pheno <- sqrt(1-pi)*rnorm(n=n, mean=0, sd=1)+start_geno[1,]*sqrt(pi/(2*q*(1-q)))
    pheno <- matrix(pheno)
    colnames(pheno)=c("pheno_trait1")
    write.table(pheno, file="sim_pheno.txt", row.names=T, col.names=T, quote=F,sep="\t")

    for(i in 1: snp)
	geno_type1[i,] <- start_geno[1,]
	
    change <- rho * n
    change <- ifelse(change%%2==0,change, change+1)

    for(i in 2: snp) {
	idx = sample(1:n, change, replace=F)
	for(j in seq(1,length(idx),2))
	    swap(geno_type1[i,idx[j]], geno_type1[i,idx[j+1]])
    }

    write.table(geno_type1, file="sim_geno1.txt", row.names=F, col.names=F, quote=F,sep="\t")

    geno_type2 <- geno_type1 - 1
    write.table(geno_type2, file="sim_geno2.txt", row.names=F, col.names=F, quote=F,sep="\t")
}
