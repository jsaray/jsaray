

rnaSilSpec = readLines('/Users/jsaray/Internship/Code/debug/silhouette_spectral-rna.txt')
rnaSilAgglo = readLines('/Users/jsaray/Internship/Code/debug/silhouette_agglom-rna.txt')

mirnaSilSpec = readLines('/Users/jsaray/Internship/Code/debug/silhouette_spectral-mirna.txt')
mirnaSilAgglo = readLines('/Users/jsaray/Internship/Code/debug/silhouette_agglom-mirna.txt')

methSilSpec = readLines('/Users/jsaray/Internship/Code/debug/silhouette_spectral-meth.txt')
methSilAgglo = readLines('/Users/jsaray/Internship/Code/debug/silhouette_agglom-meth.txt')


plot(seq(2,length(rnaSilSpec)+1),as.numeric(rnaSilSpec))
plot(seq(2,length(rnaSilAgglo)+1),as.numeric(rnaSilAgglo))

plot(seq(2,length(mirnaSilSpec)+1),as.numeric(mirnaSilSpec))
plot(seq(2,length(mirnaSilAgglo)+1),as.numeric(mirnaSilAgglo))

plot(seq(2,length(methSilSpec)+1),as.numeric(methSilSpec))
plot(seq(2,length(methSilAgglo)+1),as.numeric(methSilAgglo))


rnaElb = readLines('/Users/jsaray/Internship/Code/scriptPython/debug/elbow_kmeans_rna.txt')

mirnaElb= readLines('/Users/jsaray/Internship/Code/scriptPython/debug/elbow_kmeans_mirna.txt')

methElb = readLines('/Users/jsaray/Internship/Code/scriptPython/debug/elbow_kmeans_meth.txt')

plot(seq(2,length(rnaElb)+1),as.numeric(rnaElb))
plot(seq(2,length(mirnaElb)+1),as.numeric(mirnaElb))
plot(seq(2,length(methElb)+1),as.numeric(methElb))
