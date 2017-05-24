
library(ggfortify)


fileName = "./BatchEffect/GA.csv"
dataFrameGA = read.csv("./BatchEffect/GA.csv",sep="\t",stringsAsFactors = FALSE)
dataFrameHISEQ = read.csv("./BatchEffect/HISEQ.csv",sep="\t",stringsAsFactors = FALSE)
print(dim(dataFrameGA))
numColsGA = dim(dataFrameGA)[2]
numColsHISEQ = dim(dataFrameHISEQ)[2]

print('numCols GA is ')
print(numColsGA)
print('numCols HISEQ is ')
print(numColsHISEQ)

# Tidying Up, GA Table.
gaV1 = seq(2,763,3)
gaV2 = seq(3,763,2)
#
dataFrameGA = dataFrameGA[-gaV1]
dataFrameGA = dataFrameGA[-gaV2]
dataFrameGA = dataFrameGA[-1,]
#
# Tidying Up, HISEQ Table.
hiseqV1 = seq(2,892,3)
hiseqV2 = seq(3,892,2)
#
dataFrameHISEQ = dataFrameHISEQ[-hiseqV1]
dataFrameHISEQ = dataFrameHISEQ[-hiseqV2]
dataFrameHISEQ = dataFrameHISEQ[-1,]

dataFrameGATranspose = t(dataFrameGA[,2:ncol(dataFrameGA)])
dataFrameHISEQTranspose = t(dataFrameHISEQ[,2:ncol(dataFrameHISEQ)])

colnames(dataFrameGATranspose) = dataFrameGA[,1]
colnames(dataFrameHISEQTranspose) = dataFrameHISEQ[,1]
#
print(dim(dataFrameGA))
print(dim(dataFrameHISEQ))

ve = seq(1:705)

asNumericGA = as.numeric(dataFrameGATranspose)
matrGA = matrix(asNumericGA,nrow=254,ncol=705)
#varZeroColsGA = which(apply(matrGA, 2, var)==0)
#maskGA = !(ve %in% varZeroColsGA)
#matrGA = matrGA[,maskGA]
batchGAColumn = rep("GA",254)
matrGA = as.data.frame(matrGA)
matrGA$Batch = batchGAColumn

asNumericHISEQ = as.numeric(dataFrameHISEQTranspose)
matrHISEQ = matrix(asNumericHISEQ,nrow=297,ncol=705)
#varZeroColsHISEQ = which(apply(matrHISEQ, 2, var)==0)
#maskHISEQ = !(ve %in% varZeroColsHISEQ)
#matrHISEQ = matrHISEQ[,maskHISEQ]
batchHISEQColumn = rep("HISEQ",297)
matrHISEQ = as.data.frame(matrHISEQ)
matrHISEQ$Batch = batchHISEQColumn


PCAReady = rbind(matrGA,matrHISEQ)
batchColumn = PCAReady$Batch
varZeroColsGA = which(apply(PCAReady[,1:705], 2, var)==0)
mask = !(ve %in% varZeroColsGA)
PCAReady= PCAReady[,mask]

# At this point we have GA and HISEQ table with only reads_per_million_miRNA label

# We will PCA Plot HISEQ and GA tables in order to check if a Bias exist :
# A common method for visualizing the existence of batch effects is PCA. 
# The first two principal components are plotted with each sample colored 
# by the suspected batch, and separation of colors is taken as evidence of a batch effect.
# Source : 
# https://academic.oup.com/bioinformatics/article/29/22/2877/313226/A-new-statistic-for-identifying-batch-effects-in

#

pca = prcomp(PCAReady[,1:635],center=TRUE,scale.=TRUE)
autoplot(pca,data=PCAReady,colour='Batch')
# scaledGA = scale(dataFrameGA)
# scaledHISEQ = scale(dataFrameHISEQ)
# #
#
# #
# pcaHISEQ = prcomp(scaledHISEQ,mean=TRUE,scale.=TRUE)
# #
# print(pcaGA)
# print(pcaHISEQ)