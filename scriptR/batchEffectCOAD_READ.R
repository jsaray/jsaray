library(ggfortify)
library(mice)
library(VIM)
library(car)
coad = read.csv('/Users/jsaray/BatchEffect/COAD.miRseq_mature_RPM.csv',sep='\t',stringsAsFactors = FALSE)
reado = read.csv('/Users/jsaray/BatchEffect/READ.miRseq_mature_RPM.csv',sep='\t',stringsAsFactors = FALSE)

coadlog2 = read.csv('/Users/jsaray/BatchEffect/COAD.miRseq_mature_RPM_log2.csv',sep='\t',stringsAsFactors = FALSE)
readlog2 = read.csv('/Users/jsaray/BatchEffect/READ.miRseq_mature_RPM_log2.csv',sep='\t',stringsAsFactors = FALSE)

ve = seq(1:770)
# First PCA COAD and READ without LOG 2 ######
# 221 X 2588 DImensions of COAD

coadT = t(coad[,2:ncol(coad)]    )
readT = t(  reado[,2:ncol(reado)]   )

asNumericCoad = as.numeric(coadT)
matrCoad = matrix(asNumericCoad,nrow=221,ncol=2588)
batchCoadColumn = rep("COAD",221)
matrCoad = as.data.frame(matrCoad)
matrCoad$Batch = batchCoadColumn

asNumericRead = as.numeric(readT)
matrRead = matrix(asNumericRead,nrow=76,ncol=2588)
batchReadColumn = rep("READ",76)
matrRead = as.data.frame(matrRead)
matrRead$Batch = batchReadColumn

PCAReady = rbind(matrCoad,matrRead)
print('PCAReady Dim')
print(dim(PCAReady))
batchColumn = PCAReady$Batch
varZeroCols = which(apply(PCAReady[,1:769], 2, var)==0)

mask = !(ve %in% varZeroCols)
PCAReady= PCAReady[,mask]


print('Number of Zero Variance columns : ')
print(length(varZeroCols))

print('Disappearing all NA columns : ')
PCAReady = PCAReady[, colSums(is.na(PCAReady)) < nrow(PCAReady) - 1]
print('Purging columns with more than 5% of NaN')
v = apply(PCAReady,2,countNA)
PCAReady = PCAReady[, which( v < 11 )]

print("Dimensions after ")
print(dim(PCAReady))

print("Imputation - Dirty mean method for now, to correct later.")
column_means = colMeans(PCAReady[,1:380],na.rm=TRUE)
for (i in 1:380){
  inds = which(is.na(PCAReady[,i]))
  print('inds')
  print(inds)
  PCAReady[inds,i] <- column_means[i]
  print('pca ready aftert')
  print(PCAReady[inds,i])
}

patt = md.pattern(PCAReady)
isComplete  = any(is.na(PCAReady[,1:380]))
print('Complete Matrix')
print(isComplete)
print('PCA computing')
pca = prcomp(PCAReady[,1:380],center=TRUE,scale.=TRUE)
autoplot(pca,data=PCAReady,colour='Batch')