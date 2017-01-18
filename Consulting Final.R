####This script takes Matrix eQTL input files and creates a weight file and a predictive performance (R2) file
####The weight output file can be used to generate a .db file for use in PrediXcan with generate_sqlite_dbs.py
####by Heather E. Wheeler 20160803####

###This is script is based on Heather Wheeler's original code to run an elastic net model
###The script is modified to run SVM-KNN-PCA-NS

###The initial part of this code is provided by Dr. Wheeler from her previous study in order to create
###The data set to be used in the study.
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
args <- c('22','1','YRI') #uncomment for testing in RStudio
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables

exp.dir <- "C:/Users/Nick/Dropbox/STAT_consulting_data/"
snp.dir <- "C:/Users/Nick/Dropbox/STAT_consulting_data/"
snp.annot.dir <- "C:/Users/Nick/Dropbox/STAT_consulting_data/"
out.dir <- "C:/Users/Nick/Dropbox/STAT_consulting_data/output/"

k <- 10 ### k-fold CV
n <- 1 #number of k-fold CV replicates

##alpha = The elasticnet mixing parameter, with 0≤α≤ 1. 
#alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.

chromosome <- args[1]
alpha <- as.numeric(args[2]) #alpha to test in CV
pop <- args[3]

### alter filenames according to how you named them
exp.file <- exp.dir %&% pop %&% "_Expression.txt.gz"
exp.annot.file <- exp.dir %&% "GRCh37_hg19_ILMN_Human-6_v2_gene_annotation_for_elasticNet.txt"
snp.file <- snp.dir %&% pop %&% "_" %&% chromosome %&% ".SNP.txt.gz"
snp.annot.file <- snp.annot.dir %&% pop %&% "_" %&% chromosome %&% ".SNP.Location.txt.gz"

################################################
### Functions & Libraries
library(glmnet)
library(dplyr)
################################################

##get gene pos info
gencode <- read.table(exp.annot.file,header=TRUE)
##get snp pos info
snpcode <- read.table(snp.annot.file,header=TRUE)
##get snp allele info (needed for weight output)
allelecode <- read.table(snp.dir %&% "chr" %&% chromosome %&% "_" %&% pop %&% "_alleles.txt.gz")
colnames(allelecode) <- c("CHR","POS","SNP","refAllele","effectAllele")
rownames(allelecode) <- allelecode$POS #name by position b/c we found duplicate rsids

##read exp and chr gt dosages
exp <- read.table(exp.file, header=TRUE)
gt <- read.table(snp.file,header=TRUE)

##join pos info 
popgt <- left_join(snpcode,gt,by=c("snp"="id"))
popgt <- popgt[duplicated(popgt$snp)==FALSE,] #remove duplicated rsids with incorrect pos
popgt <- popgt[duplicated(popgt$pos)==FALSE,] #remove duplicated pos 
rownames(popgt) <- popgt[,3] #name by position b/c we found duplicate rsids
popgt <- popgt[popgt[,3] %in% allelecode$POS,] #only keep SNPs in allelecode file (removes esv SNPs)
##join gene info
popexp <- left_join(gencode,exp,by=c("geneid"="id"))

popsamplelist <- colnames(exp)[-1]

#pull gene info & expression from pop of interest
popexp <- dplyr::filter(popexp,chrom==chromosome)
explist <- as.character(popexp$geneid)

set.seed(42)
groupid <- sample(1:10,length(popsamplelist),replace=TRUE) ##need to use same folds to compare alphas

resultsarray <- array(0,c(length(explist),8))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")
dimnames(resultsarray)[[2]] <- resultscol
workingbest <- out.dir %&% "HapMap3_" %&% pop %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_chr" %&% chromosome %&% "_" %&% date %&% ".txt"
write(resultscol,file=workingbest,ncolumns=8,sep="\t")

weightcol = c("gene", "rsid", "ref", "alt", "beta", "alpha") #col headers for use with generate_sqlite_dbs.py
workingweight <- out.dir %&% "HapMap3_" %&% pop %&% "_elasticNet_alpha" %&% alpha %&% "_weights_chr" %&% chromosome %&% "_" %&% date %&% ".txt"
write(weightcol,file=workingweight,ncolumns=6,sep="\t")

R2.01<-matrix()
A.R2.01<-matrix()
R2.11<-matrix()
A.R2.11<-matrix()
R2.model1<-matrix()

library(e1071)
library(kknn)
library(splines)

for(i in 1:length(explist)){
 cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  start <- popexp$s1[i] - 1e6 ### 1Mb gene lower bound for cis-eQTLS
  end <- popexp$s2[i] + 1e6 ###  1Mb gene upper bound for cis-eQTLs
  cisgenos <- subset(popgt,popgt[,3]>=start & popgt[,3]<=end) ### pull cis-SNP genotypes
  rownames(cisgenos) <- cisgenos$pos #carry positions along
  cismat <- as.matrix(cisgenos[,4:dim(cisgenos)[2]]) #get dosages only in matrix format for glmnet
  cismat <- t(cismat) #transpose to match previous code
  expmat <- as.matrix(popexp[,9:dim(popexp)[2]]) #make exp only in matrix format for glmnet
  expmat <- t(expmat) #transpose to match previous code
  colnames(expmat) <- popexp$geneid #carry gene IDs along
  if(is.null(dim(cismat))){
    #if(is.null(dim(cismat)) | gene=="ILMN_1740816"){ #special case for GIH alpha=0 to skip Error in predmat[which, seq(nlami)] = preds : replacement has length zero
    bestbetas <- data.frame() ###effectively skips genes with 0 cis-SNPs
  }else{
    minorsnps <- subset(colMeans(cismat), colMeans(cismat,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
    minorsnps <- names(minorsnps)
    cismat <- cismat[,minorsnps]
    if(length(minorsnps) < 2){###effectively skips genes with <2 cis-SNPs
      bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{
      exppheno <- expmat[,gene] ### pull expression data for gene
      exppheno <- scale(exppheno, center=T, scale=T)  ###scale to compare across genes
      exppheno[is.na(exppheno)] <- 0

#Code that I was responsible for writing with the creation of the data set 
      colnames(exppheno) <- "exp"
      pe <- cbind(exppheno, cismat)
      pe <- as.data.frame(pe)
      pe.best.k <- cbind(exppheno,cismat)
      pe.best.k <- as.data.frame(pe.best.k)
      pe.best.k$exp <- as.factor(ifelse(pe.best.k$exp > 0, "1", "0"))
     
      train.k<-train.kknn(exp~.,pe.best.k,ks=c(1:15),distance=1)
      best.k<-train.k$best.parameters$k
   
      
     #Runs the model with the best k
      
      svm.data<-pe[,-max(dim(pe.best.k)[2])]
      model.svm<-svm(exp~.,data=svm.data,scale=FALSE)
      final.train<-pe.best.k[model.svm$index,-max(dim(pe.best.k)[2])]
      final.test<-pe.best.k[-max(dim(pe.best.k)[2])]
      model.knn<-kknn(exp~.,final.train,final.test,k=best.k,kernel="rectangular",dist=1)
      
      fv <- model.knn$fitted.values
      gp1 <- which(fv==0)
      dat0 <- as.data.frame(pe[gp1,])
      dat00<-as.data.frame(pe[gp1,-1])
      dat1 <- as.data.frame(pe[-gp1,])
      dat01<-as.data.frame(pe[-gp1,-1])
      
      pc0 <- prcomp(t(dat00))
      pcom0 <- pc0$rotation
      z<-summary(pc0)
      l<-z$importance[2,]
      m<-which(l>=.1)
      co0<-max(m)
      
      st0 <- lm(paste("dat0$exp ~", paste(paste0("ns(pcom0[,", 1:co0, "],3)"), collapse = "+")))
      mod.sum0 <-summary(st0)
      
      pc1 <- prcomp(t(dat01))
      pcom1 <- pc1$rotation
      z1<-summary(pc1)
      l1<-z1$importance[2,]
      m1<-which(l1>=.1)
      co1<-max(m1)
      
      st1 <- lm(paste("dat1$exp ~", paste(paste0("ns(pcom1[,", 1:co1, "],3)"), collapse = "+")))
      
      mod.sum1 <- summary(st1)
      
      ybar0<-mean(dat0$exp)
      ybar1<-mean(dat1$exp)
      SSE<-sum(sum((st0$fitted.values-dat0$exp)**2),sum((st1$fitted.values-dat1$exp)**2))
      SST<-sum(sum((dat0$exp-ybar0)**2),sum((dat1$exp-ybar1)**2))
      R2.f<-1-(SSE/SST)
      }
  }
  R2.01[i] <- mod.sum0$r.squared
  A.R2.01[i] <- mod.sum0$adj.r.squared
  
  R2.11[i] <- mod.sum1$r.squared
  A.R2.11[i] <- mod.sum1$adj.r.squared
  
  R2.model1[i]<-R2.f
}
