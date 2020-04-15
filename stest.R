# check if package Kendall has been installed
if(!require(Kendall)) install.packages("Kendall")
library(Kendall)
    
stest <- function(geno1, geno2, pheno, n_pca) {
    if(is.null(geno1))	stop("Error: provide genotype encoded as E1")
    if(is.null(geno2))	stop("Error: provide genotype encoded as E2")
    if(is.null(pheno))	stop("Error: provide phenotype data")
    if(is.null(n_pca)) stop("Error: enter number of principal components")

    geno_type1 <- read.table(geno1)
    geno_type1 <- as.matrix(geno_type1)

    geno_type2 <- read.table(geno2)
    geno_type2 <- as.matrix(geno_type2)

    print("genotype encoding completed.")

    # loading phenotype data 
    pheno<-read.table(pheno,header=T,stringsAsFactor=F)
    pheno<-as.matrix(pheno[,-1])
    if(ncol(pheno)< 2)
         colnames(pheno) <- c("pheno_trait1")
    
    # performing pca 
    pc_type1 <- prcomp(geno_type1, center=TRUE, scale=TRUE) 
    pc_type2 <- prcomp(geno_type2, center=TRUE, scale=TRUE) 
	            
    n_pca <- as.numeric(n_pca)    
    pca_vec_type1 <- pc_type1$rotation # eigenvectors 
    pca_vec_type2 <- pc_type2$rotation # eigenvectors 

    print("PCA completed.")

    # genotype adjustment 
    geno_type1_2 <- geno_type1 %*% pca_vec_type1[,1:n_pca] 
    geno_type1_3 <- geno_type1_2 %*% t(pca_vec_type1[,1:n_pca])
    adj_geno_type1 <- geno_type1 - geno_type1_3
    
    geno_type2_2 <- geno_type2 %*% pca_vec_type2[,1:n_pca] 
    geno_type2_3 <- geno_type2_2 %*% t(pca_vec_type2[,1:n_pca])
    adj_geno_type2 <- geno_type2 - geno_type2_3
    
    print("genotypes have been adjusted")
    
    kendall_skewness <- c(1:ncol(pheno))
    kendall_kurtosis <- c(1:ncol(pheno))
    pearson_skewness <- c(1:ncol(pheno))
    pearson_kurtosis <- c(1:ncol(pheno))

    for (i in 1:ncol(pheno)) {
        
        # Kendall's test with E1
        kendall_pvalue_type1 <- c(1:nrow(geno_type1))
        for(n in 1:nrow(geno_type1))
             kendall_pvalue_type1[n]= Kendall(adj_geno_type1[n,], pheno[,i])$sl

        # Kendall's test with E2
        kendall_pvalue_type2 <- c(1:nrow(geno_type2))
        for(n in 1:nrow(geno_type2))
             kendall_pvalue_type2[n]= Kendall(adj_geno_type2[n,], pheno[,i])$sl
         
        # Pearson's test with E1	
        pearson_pvalue_type1 <-c(1:nrow(geno_type1))
        for(n in 1:nrow(geno_type1))
	     pearson_pvalue_type1[n]=(cor.test(adj_geno_type1[n,], pheno[,i], method="pearson"))$p.value
     
        # Pearson's test with E2	
        pearson_pvalue_type2 <-c(1:nrow(geno_type2))
        for(n in 1:nrow(geno_type2))
             pearson_pvalue_type2[n]=(cor.test(adj_geno_type2[n,], pheno[,i], method="pearson"))$p.value
    
        # get difference in p-values (delta p-value) between E1 and E2	 
        kendall_diff <- kendall_pvalue_type1 - kendall_pvalue_type2
        pearson_diff <- pearson_pvalue_type1 - pearson_pvalue_type2		  
    
        # plot histogram of delta p-value
        pheno_name<- colnames(pheno)
        file_name=paste(pheno_name[i],".jpg")
        file_name=gsub(" ","", file_name)
        jpeg(filename=file_name, width=800, height=400)
        par(mfrow=c(1,2))	 
        hist(kendall_diff, freq=F, main="Kendall", xlab="delta p", ylab="Relative frequency", breaks=20)
        hist(pearson_diff, freq=F, main="Pearson", xlab="delta p", ylab="Relative frequency", breaks=20)
        dev.off()
	 
        # estimate skewness and kurtosis
        n <- length(kendall_diff)

        kendall_m3 <- sum(kendall_diff^3)/n
        kendall_m2 <- sum(kendall_diff^2)/n
        kendall_m4 <- sum(kendall_diff^4)/n

        kendall_skewness[i] <- (sqrt(n*(n-1))*kendall_m3)/((n-2)*(kendall_m2)^(3/2))
        kendall_kurtosis[i] <- ((n-1)*(n+1)*kendall_m4)/((n-2)*(n-3)*kendall_m2^2) - (3*(n-1)^2)/((n-2)*(n-3))
    
        pearson_m3 <- sum(pearson_diff^3)/n
        pearson_m2 <- sum(pearson_diff^2)/n
        pearson_m4 <- sum(pearson_diff^4)/n
       
        pearson_skewness[i] <- (sqrt(n*(n-1))*pearson_m3)/((n-2)*(pearson_m2)^(3/2))
        pearson_kurtosis[i] <- ((n-1)*(n+1)*pearson_m4)/((n-2)*(n-3)*pearson_m2^2) - (3*(n-1)^2)/((n-2)*(n-3))

   }

    skew<- cbind(kendall_skewness, pearson_skewness)
    kurt<- cbind(kendall_kurtosis, pearson_kurtosis)


    jpeg(filename="boxplots.jpg", width=800, height=400)
    par(mfrow=c(1,2))	 
    boxplot(skew,col=c("orange","white"), border=c("black","black"), main="Skewness", names=c("Kendall","Pearson"))
    boxplot(kurt,col=c("orange","white"), border=c("black","black"), main="Kurtosis", names=c("Kendall","Pearson"))
    dev.off()
    
    print("stest completed")
}
# end of stest.R

