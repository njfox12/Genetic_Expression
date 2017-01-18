###This code takes matrix eQTL input files and runs K-Nearest Neighbors (KNN), 
###Principal Component Analysis (PCA), and natural cubic splines (NS). 
###By Stephanie Oliva
###12/17/2016

library(kknn)
library(splines)
library(bigpca)
library(dplyr)

###############################################
##This section of the code was provided by Dr. Wheeler.
##It has been edited to include only what's needed for the KNN-PCA-NS model.

date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
args <- c('22','1','YRI') #uncomment for testing in RStudio
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables

exp.dir <- "~/Dropbox/STAT_consulting_data/"
snp.dir <- "~/Dropbox/STAT_consulting_data/"
snp.annot.dir <- "~/Dropbox/STAT_consulting_data/"
out.dir <- "~/Dropbox/"

chromosome <- args[1]
alpha <- as.numeric(args[2]) 
pop <- args[3]

### alter filenames according to how you named them
exp.file <- exp.dir %&% pop %&% "_Expression.txt.gz"
exp.annot.file <- exp.dir %&% "GRCh37_hg19_ILMN_Human-6_v2_gene_annotation_for_elasticNet.txt"
snp.file <- snp.dir %&% pop %&% "_" %&% chromosome %&% ".SNP.txt.gz"
snp.annot.file <- snp.annot.dir %&% pop %&% "_" %&% chromosome %&% ".SNP.Location.txt.gz"


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

##End of this section of Dr. Wheeler's code.
###################################################

##vectors to be updated for subset categorized as 0
R2.0 <- NULL #R2 
A.R2.0 <- NULL #adj R2
comp0 <- NULL #number of components 
n.0 <- NULL #number of obs
p.0 <- NULL #total predictors in NS model

##vectors to be updated for subset categorized as 1
R2.1 <- NULL #R2
A.R2.1 <- NULL #adj R2
comp1 <- NULL #number of components
n.1 <- NULL #number of obs
p.1 <- NULL #total predictors in NS model

k_fin <- NULL #best K per gene
g <- NULL #gene name
r2.fin <- NULL #combined R2 for each gene

#################################################
##First part of this loop written by Dr. Wheeler
################################################

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
    }
  }
##########################################
##End of Dr. Wheeler's code
##########################################
colnames(exppheno) <- "expr" #naming response
pe <- cbind(exppheno, cismat) #binding response and predictors for KNN
cvdat <- as.data.frame(pe)
cvdat$expr <- ifelse(cvdat$expr < 0, 0, 1) #categorizing response as 0,1

########################################
############ KNN #######################

cv <- train.kknn(as.factor(expr) ~., data = cvdat, ks=c(2:15), 
                 distance=2, kernel = "rectangular") #LOOCV for knn (much faster)
k <- unlist(cv$best.parameters[2]) 
k_fin[i] <-  k #saving best K
kn1 <- kknn(as.factor(expr)~., cvdat, cvdat, k=k, 
            dist=2, kernel = "rectangular") #fitting knn with best K
#####dist=2 is Euclidean, dist=1 was run for Manhattan

fv <- kn1$fitted.values #saving fitted values

###subsetting data by cat 0, 1
gp1 <- which(fv==0) 
dat0 <- as.data.frame(pe[gp1,]) 
dat1 <- as.data.frame(pe[-gp1,]) 

##################################
###########PCA cat 0##############

pc0 <- prcomp(t(dat0[,-1])) #running PCA
pcom0 <- pc0$rotation #principal components
co0 <- quick.elbow(pc0$sdev^2) #finding principal components that explain the 
#most variability

##################################
###########NS Model cat 0#########

st0 <- lm(paste("dat0$expr ~", paste(paste0("ns(pcom0[,", 1:co0, "],3)"), 
                                     collapse = "+"))) #fitting NS model with 3 df
mod.sum0 <-summary(st0) 
R2.0[i] <- mod.sum0$r.squared #retrieving R2 values
A.R2.0[i] <- mod.sum0$adj.r.squared #adj R2 values
comp0[i] <- quick.elbow(pc0$sdev^2) #number of principal components

##################################
###########PCA cat 1##############

pc1 <- prcomp(t(dat1[,-1])) #PCA on subset with cat 1
pcom1 <- pc1$rotation #principal components
co1 <- quick.elbow(pc1$sdev^2) #finding principal components that explain the
#highest variability

##################################
###########NS Model cat 1#########

st1 <- lm(paste("dat1$expr ~", paste(paste0("ns(pcom1[,", 1:co1, "],3)"),
                                     collapse = "+"))) #NS model on pc

mod.sum1 <- summary(st1)
R2.1[i] <- mod.sum1$r.squared #R2 values
A.R2.1[i] <- mod.sum1$adj.r.squared #adj R2 values
comp1[i] <- quick.elbow(pc1$sdev^2) #number of pc

##################################
#####Calculating combined R2######

###predicted values
pred0 <- predict(st0) 
pred1 <- predict(st1) 

###mean response 
ybar.0 <- mean(dat0$expr) 
ybar.1 <- mean(dat1$expr)
ybar.fin <- c(rep(ybar.0, nrow(dat0)), rep(ybar.1, nrow(dat1))) #vector of means
dat.p <- rbind(dat0, dat1) #combine data for response
pred.p <- c(pred0, pred1) #combine predictions

###R2
r2.fin[i] <- 1- (sum((dat.p$expr - pred.p)^2) / sum((dat.p$exp-ybar.fin)^2)) 

##################################
#####Saving values of interest####

g[i] <- gene #keeping gene name

###n obs
n.0[i] <- nrow(dat0) 
n.1[i] <- nrow(dat1) 

###NS predictors
p.0[i] <- co0 *3 
p.1[i] <- co1 *3 
}


fin <- cbind(g,k_fin, comp0, n.0, p.0, R2.0, A.R2.0, comp1, n.1,
             p.1, R2.1, A.R2.1, r2.fin) #binding data
write.csv(fin, "~/Dropbox/fin_res_k2_15_dist2_euc_finalized.csv") #saving file
