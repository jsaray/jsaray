
# normalizedRNA = read.csv("/Users/jsaray/ApplicationsJose/stddata__2016_01_28/COADREAD/20160128/ga_RSEM_genes_normalized__data/ga_RSEM_genes_normalized__data.csv",sep="\t",stringsAsFactors = FALSE)
# normalizedRNA = t(normalizedRNA)
# print(dim(normalizedRNA))
# colnames(normalizedRNA) = normalizedRNA[1,]
# print(colnames(normalizedRNA))
# geneA1BG = normalizedRNA[,c("AANAT|15")]
# geneA1BG = as.numeric( geneA1BG[2:length(geneA1BG)] )
# hist(geneA1BG)

df = read.csv("/Users/jsaray/DataFrameCreation/miRseq_mature_RPM.csv",sep="\t",stringsAsFactors = FALSE)
dflog = read.csv("/Users/jsaray/DataFrameCreation/miRseq_mature_RPM_log2.csv",sep="\t",stringsAsFactors = FALSE)

dflog = t(dflog )
print(dim(dflog ))
colnames(dflog ) = dflog [1,]
normalities = c()
print(dim(dflog))
for(i in 2:2588){
print(i)
  if (  length(which(dflog[,i] != 'NA')) >= 4)
  {
    numbers = as.numeric ( dflog[,i] )
    normalities = append(normalities,shapiro.test(numbers)$p.value)
  }
}

co = which(normalities > 0.05)
print(length(co))
#logDf = log( as.numeric(dflog [ 1:298,1:2588 ]     ) )
#logDf = matrix(dflog ,nrow=298,ncol=2588)
#print (logDf[,1])