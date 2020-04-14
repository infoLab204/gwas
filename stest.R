if(!require(Kendall)) install.packages("Kendall")
library(Kendall)
    
stest <- function(geno1,geno2,pheno, pca_num) {
    if(is.null(pheno))	stop("phenotypes must exist.")
    if(is.null(geno1))	stop("First encoded genotype")
    if(is.null(geno2))	stop("Second encoded genotype")
    if(is.null(pca_num))	stop("please,pca_num")

    geno_type1 <- read.table(geno1)
    geno_type1 <- as.matrix(geno_type1)

    geno_type2 <- read.table(geno2)
    geno_type2 <- as.matrix(geno_type2)

    print("genoenotype encoding complete.")

    # phenotype data loading
    pheno<-read.table(pheno,header=T,stringsAsFactor=F)
    pheno<-as.matrix(pheno[,-1])
    
    # pca 
    pc_type1 <- prcomp(geno_type1, center=TRUE, scale=TRUE) 
    pc_type2 <- prcomp(geno_type2, center=TRUE, scale=TRUE) 
	            
    pca_num <- as.numeric(pca_num)    
    pca_vec_type1 <- pc_type1$rotation # eigenvectors 
    pca_vec_type2 <- pc_type2$rotation # eigenvectors 

    print("PCA complete.")

    #adjusted genotype 
    geno_type1_2 <- geno_type1 %*% pca_vec_type1[,1:pca_num] 
    geno_type1_3 <- geno_type1_2 %*% t(pca_vec_type1[,1:pca_num])
    adj_geno_type1 <- geno_type1 - geno_type1_3
    
    geno_type2_2 <- geno_type2 %*% pca_vec_type2[,1:pca_num] 
    geno_type2_3 <- geno_type2_2 %*% t(pca_vec_type2[,1:pca_num])
    adj_geno_type2 <- geno_type2 - geno_type2_3
    
    print("Adjusted genoenotype complete.")
    
    kendall_skewness <- c(1:ncol(pheno))
    kendall_kurtosis <- c(1:ncol(pheno))
    pearson_skewness <- c(1:ncol(pheno))
    pearson_kurtosis <- c(1:ncol(pheno))

#for (i in 1:1) {
for (i in 1:ncol(pheno)) {
        
        # kendall test type1
        kendall_pvalue_type1 <- c(1:nrow(geno_type1))
        for(n in 1:nrow(geno_type1))
             kendall_pvalue_type1[n]= Kendall(adj_geno_type1[n,], pheno[,i])$sl

        # kendall test type2
        kendall_pvalue_type2 <- c(1:nrow(geno_type2))
        for(n in 1:nrow(geno_type2))
             kendall_pvalue_type2[n]= Kendall(adj_geno_type2[n,], pheno[,i])$sl
         
        # pearson test type1	
        pearson_pvalue_type1 <-c(1:nrow(geno_type1))
        for(n in 1:nrow(geno_type1))
	     pearson_pvalue_type1[n]=(cor.test(adj_geno_type1[n,], pheno[,i], method="pearson"))$p.value
     
        # pearson test type2	
        pearson_pvalue_type2 <-c(1:nrow(geno_type2))
        for(n in 1:nrow(geno_type2))
             pearson_pvalue_type2[n]=(cor.test(adj_geno_type2[n,], pheno[,i], method="pearson"))$p.value
    
        kendall_diff <- kendall_pvalue_type1 - kendall_pvalue_type2
        pearson_diff <- pearson_pvalue_type1 - pearson_pvalue_type2		  
    
        # histogram of pvalue's difference
        pheno_name<- colnames(pheno)
        file_name=paste(pheno_name[i],".jpg")
        file_name=gsub(" ","", file_name)
        jpeg(filename=file_name, width=800, height=400)
        par(mfrow=c(1,2))	 
        hist(kendall_diff, freq=F, main="Kendall", xlab="delta p", ylab="Relative frequency", breaks=20)
        hist(pearson_diff, freq=F, main="Pearson", xlab="delta p", ylab="Relative frequency", breaks=20)
        dev.off()
	 
        n <- length(kendall_diff)

        kendall_m3 <- sum(kendall_diff^3)/n
        kendall_m2 <- sum(kendall_diff^2)/n
        kendall_m4 <- sum(kendall_diff^4)/n

        kendall_skewness[i] <- (sqrt(n*(n-1))/(n-2))*(kendall_m3/((kendall_m2)^(3/2)))
        kendall_kurtosis[i] <- ((n-1)/((n-2)*(n-3)))*{(n+1)*(kendall_m4/(kendall_m2)^2)}-3
    
        pearson_m3 <- sum(pearson_diff^3)/n
        pearson_m2 <- sum(pearson_diff^2)/n
        pearson_m4 <- sum(pearson_diff^4)/n
        
        pearson_skewness[i] <- (sqrt(n*(n-1))/(n-2))*(pearson_m3/((pearson_m2)^(3/2)))
        pearson_kurtosis[i] <- ((n-1)/((n-2)*(n-3)))*{(n+1)*(pearson_m4/(pearson_m2)^2)}-3

   }

    skew<- cbind(kendall_skewness, pearson_skewness)
    kurt<- cbind(kendall_kurtosis, pearson_kurtosis)
    print(skew)
    print(kurt)

    jpeg(filename="boxplots.jpg", width=800, height=400)
    par(mfrow=c(1,2))	 
    boxplot(skew,col=c("orange","white"), border=c("black","black"), main="Skewness", names=c("Kendall","Pearson"))
    boxplot(kurt,col=c("orange","white"), border=c("black","black"), main="Kurtosis", names=c("Kendall","Pearson"))
    dev.off()
    
    print("GWAS complete.")
}

