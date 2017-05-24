install.packages('iCluster')
library('iCluster')

rna = read.csv('/Users/jsaray/DataFrameCreation/readyToMixRNA.csv')
mirna = read.csv('/Users/jsaray/DataFrameCreation/readyToMixMIRNA.csv')
meth = read.csv('/Users/jsaray/DataFrameCreation/readyToMixMETH.csv')

rna = rna[,2:ncol(rna)]
mirna = mirna[,2:ncol(mirna)]
meth  = meth[,2:ncol(meth)]

toCluster = list(as.matrix(rna),as.matrix(mirna),as.matrix(meth))

fit = iCluster(toCluster, k=4,lambda=c(0.2,0.2,0.2))
#plotiCluster(fit=fit,l)
