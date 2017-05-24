

install.packages('clusteval')
install.packages('rgl')
library(clusteval)
library(kernlab)
library(car)
library(Rcmdr)
library(clusteval)
library(rgl)
library(omicade4)

# This script will perform the following :
# * Takes three files containing the same number of patients and different omics information.
# * Perform value Imputation.
# Function transformations ?.
# Apply PCA and 3 different clustering algorithms and plot to visualize
# 

applyClusteringAlgo <- function(df,nClusters,nuRows,nuCols){
  return(specc(as.matrix(df,nrow=nuRows,ncol=nuCols),centers=nClusters))
}

topNCorrelationsPerComponent <- function(matrixRot,topn){
  for (i in 1:ncol(matrixRot)){
    col = matrixRot[,i]
    filtMask = which(col >= 0.4 || col <= -0.4)
    if (length(filtMask) == 0){
      print('NOTHING')
    }else{
      print(      paste( 'SOMETHING -',i,toString(filtMask),' values ',toString(col[filtMask]) )    )
    }
  }
}

pcaAndGraph <- function(dataFrame,titleGraph,nClusters,nuRows,nuCols,plots=FALSE){
  
  print(paste('PCA and Graph for ',titleGraph))
  pcaObject = prcomp(dataFrame,center=TRUE,scale=TRUE,retx=TRUE)

  pcaResults = summary(pcaObject)
  vectorOfVariances = pcaResults[["importance"]][2,]
  write(vectorOfVariances,file=paste("/Users/jsaray/Internship/Code/debugR/compVariances_", titleGraph,'.txt'),ncolumns=length(vectorOfVariances),sep=',')
  write.table(pcaObject$rotation,file=paste("/Users/jsaray/Internship/Code/debugR/rotation_", titleGraph,'.csv'),sep=',')
  resultClustering = applyClusteringAlgo(dataFrame,nClusters,nuRows,nuCols)
  eigenMatrix = pcaObject$rotation
  print(" Dimensions eigenRna ")
  dim(eigenMatrix)
  print(" Dimensions Rna ")
  dim(dataFrame)
  
  transfo = pcaObject$x
  
  colors <- c("#999999", "#E69F00", "#56B4E9","#00FF00","#FF0000","#000000")
  colorsCluster = colors[resultClustering@.Data]
  if(plots){
    open3d()
    scatter3d(transfo[,1],transfo[,2],transfo[,3],point.col=colorsCluster,surface=FALSE,xlab=titleGraph)
  }
  list(data=resultClustering@.Data,varVector=vectorOfVariances,cor=pcaObject$rotation)
}

returnMaxMin <- function(matr,pc){
  ret = paste( "(", max(matr[,pc]) ,  "," ,  min(matr[,pc])   , ")"    )
  ret
}

vectorFromList <- function(list){
  vtoret = c()
  for (i in seq(1,length(list))){
    vtoret = append(vtoret,list[[i]])
  }
  vtoret
}

jaccardMatrix <- function(clus1,clus2,tcgaCodes){
  
  distinctClus1 = unique(clus1)
  distinctClus2= unique(clus2)
  jacMatr = c()
for (i in sort(distinctClus1)){
    
    maskClus1 = which(clus1 == i)
    set1 = tcgaCodes[maskClus1]
    for (j in sort(distinctClus2)){
      maskClus2 = which(clus2 == j)
      set2 = tcgaCodes[maskClus2]
      jacCoef = length(  intersect(set1,set2)  ) /  length(  union(set1,set2)   )
      print( paste(jacCoef,'comparing sets ' , i , ' and ' , j)    )
      jacMatr = append(jacMatr,jacCoef)
    }
    
  }
  print(jacMatr)
  matr = matrix(jacMatr,nrow=length(distinctClus1),ncol=length(distinctClus2),byrow=TRUE)
  return(matr )
}

rna= read.csv('/Users/jsaray/DataFrameCreation/readyToMixRNA.csv',sep=',')
mirna = read.csv('/Users/jsaray/DataFrameCreation/readyToMixMIRNA.csv',sep=',')
meth = read.csv('/Users/jsaray/DataFrameCreation/readyToMixMETH.csv',sep=',')

tcgaCodes = rna[,1]
# Uncomment to take the Jaccard Index = Card(A inter B) / Card ( A Union B)
# rna = read.csv('/Users/jsaray/DataFrameCreation/cleanedRna.csv',sep=',')
# mirna = read.csv('/Users/jsaray/DataFrameCreation/cleanedMiRna.csv',sep=',')
# meth = read.csv('/Users/jsaray/DataFrameCreation/cleanedMeth.csv',sep=',')

print('dim rna')
print(dim(rna))
print('dim mirna')
print(dim(mirna))
print('dim meth')
print(dim(meth))

rna = rna[1:nrow(rna),2:ncol(rna)]
diagoRNA = diag(cov(rna))
top500VarianceRNA = order(diagoRNA)[1:500]
rna = rna[,top500VarianceRNA]
print('rna new dimensions')
print(dim(rna))
rna = scale(rna)
#rna = rna[order(patients),]

mirna = mirna[1:nrow(mirna),2:ncol(mirna)]
diagoMirna = diag(cov(mirna))
top500VarianceMirna = order(diagoMirna)[1:100]
mirna = mirna[,top500VarianceMirna]
print('mirna new dimensions')
print(dim(mirna))
mirna = scale(mirna)

meth = meth[1:nrow(meth),2:ncol(meth)]

diagoMeth = diag(cov(meth))
top500VarianceMeth = order(diagoMeth)[1:500]
meth = meth[,top500VarianceMeth]
print('meth new dimensions')
print(dim(meth))
meth = scale(meth)


logrna = sign(rna)*log(abs(rna) + 1)
logmirna = sign(mirna)*log( abs(mirna) + 1)
logmeth = sign(meth)*log( abs(meth)  + 1)

cluster_rna = pcaAndGraph(rna,'Rna sequencing',6,577,500)
cluster_mirna = pcaAndGraph(mirna,'Mirna sequencing',5,288,500)
cluster_meth = pcaAndGraph(meth,'Methylation',4,441,500)

seqClRna = seq(1,length(cluster_rna$varVector))
seqClMirna= seq(1,length(cluster_mirna$varVector))
seqClMeth= seq(1,length(cluster_meth$varVector))

barplot(vectorFromList(cluster_rna$varVector[1:50]*100),main='RNA Variance Ratio')
barplot(vectorFromList(cluster_mirna$varVector[1:50]*100),main='miRNA Variance Ratio')
barplot(vectorFromList(cluster_meth$varVector[1:50]*100),main='Meth Variance Ratio')

print(paste('Max and Min correlation with PC1 in RNA ',returnMaxMin(cluster_rna$cor,1)))
print(paste('Max and Min correlation with PCA1 in miRNA ',returnMaxMin(cluster_mirna$cor,1)))
print(paste('Max and Min correlation with PCA1 in methylation ',returnMaxMin(cluster_mirna$cor,1)))

print(paste('Max and Min correlation with PC2 in RNA ',returnMaxMin(cluster_rna$cor,2)))
print(paste('Max and Min correlation with PCA2 in miRNA ',returnMaxMin(cluster_mirna$cor,2)))
print(paste('Max and Min correlation with PCA2 in methylation ',returnMaxMin(cluster_mirna$cor,2)))

print(paste('Max and Min correlation with PC3 in RNA ',returnMaxMin(cluster_rna$cor,3)))
print(paste('Max and Min correlation with PCA3 in miRNA ',returnMaxMin(cluster_mirna$cor,3)))
print(paste('Max and Min correlation with PCA3 in methylation ',returnMaxMin(cluster_mirna$cor,3)))

# print('============ LOG ========================')
# cluster_logrna = pcaAndGraph(logrna,'500MaxLog-Rna sequencing',6,577,500)
# cluster_logmirna = pcaAndGraph(logmirna,'500MaxLog Mirna sequencing',5,288,500)
# cluster_logmeth = pcaAndGraph(logmeth,'500MaxLog Methylation',4,441,500)
# 
# seqClRna = seq(1,length(cluster_logrna$varVector))
# seqClMirna= seq(1,length(cluster_logmirna$varVector))
# seqClMeth= seq(1,length(cluster_logmeth$varVector))
# 
# barplot(vectorFromList(cluster_logrna$varVector[1:50]*100),main='Log RNA Variance Ratio')
# barplot(vectorFromList(cluster_logmirna$varVector[1:50]*100),main='Log miRNA Variance Ratio')
# barplot(vectorFromList(cluster_logmeth$varVector[1:50]*100),main='Log Meth Variance Ratio')
# print(" ========== END LOG +++++++++++++++++++++++ ")
sim_rna_mirna = cluster_similarity(cluster_rna$data,cluster_mirna$data,similarity="jaccard")
sim_rna_meth = cluster_similarity(cluster_rna$data,cluster_meth$data,similarity="jaccard")
sim_mirna_meth = cluster_similarity(cluster_mirna$data,cluster_meth$data,similarity="jaccard")

print('Sim. RNA MIRNA ')
print(sim_rna_mirna)
print('Sim. RNA METH')
print(sim_rna_meth)
print('Sim. MIRNA METH ')
print(sim_mirna_meth)

print('Sim. LOG RNA LOG MIRNA ')
print(sim_logrna_logmirna)
print('Sim. LOG RNA LOG METH')
print(sim_logrna_logmeth)
print('Sim. LOG MIRNA LOG METH ')
print(sim_logmirna_logmeth)

print(' Computing Jaccard Matrix RNA ')
jac1 = jaccardMatrix(cluster_rna$data,cluster_mirna$data,tcgaCodes)
print(jac1)

print(' Computing Jaccard Matrix MIRNA')
jac2 = jaccardMatrix(cluster_rna$data,cluster_meth$data,tcgaCodes)
print(jac2)

print(' Computing Jaccard Matrix METH ')
jac3 = jaccardMatrix(cluster_mirna$data,cluster_meth$data,tcgaCodes)
print(jac3)

print( ' JACARD LOG  MATRIXES ')

print(' Computing Jaccard Matrix LOG RNA ')
jac1 = jaccardMatrix(cluster_logrna$data,cluster_logmirna$data,tcgaCodes)
print(jac1)

print(' Computing Jaccard Matrix LOG MIRNA')
jac2 = jaccardMatrix(cluster_logrna,cluster_logmeth,tcgaCodes)
print(jac2)

print(' Computing Jaccard Matrix LOG  METH ')
jac3 = jaccardMatrix(cluster_logmirna,cluster_logmeth,tcgaCodes)
print(jac3)

listOmicade = c( t(rna),t(mirna),t(meth) )
mcoin = mcia(listOmicade,cia.nf = 2, cia.scan = TRUE, nsc = TRUE, svd = TRUE)





