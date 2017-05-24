
install.packages("SNFtool")
library(SNFtool)
## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 10; 	# Number of Iterations, usually (10~20)

## Data1 is of size n x d_1, where n is the number of patients, d_1 is the number of genes, e.g.
## Data2 is of size n x d_2, where n is the number of patients, d_2 is the number of methylation, e.g.
rna = read.csv('/Users/jsaray/DataFrameCreation/readyToMixRNA.csv')
mirna = read.csv('/Users/jsaray/DataFrameCreation/readyToMixMIRNA.csv')
meth = read.csv('/Users/jsaray/DataFrameCreation/readyToMixMETH.csv')

rna = rna[,2:ncol(rna)]
mirna = mirna[,2:ncol(mirna)]
meth  = meth[,2:ncol(meth)]
truelabel = c(matrix(1,273,1),matrix(2,273,1),matrix(3,273,1)); ##the ground truth of the simulated data;

rna = standardNormalization(rna)
mirna = standardNormalization(mirna)
meth = standardNormalization(meth)

## Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows; if the data is discrete, we recommend the users to use ""
Dist1 = dist2(as.matrix(rna),as.matrix(rna));
Dist2 = dist2(as.matrix(mirna),as.matrix(mirna));
Dist3 = dist2(as.matrix(meth),as.matrix(meth));

## next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)
W3 = affinityMatrix(Dist3, K, alpha)
## These similarity graphs have complementary information about clusters.
displayClusters(W1,truelabel[1:273]);
displayClusters(W2,truelabel[274:546]);
displayClusters(W3,truelabel[547:819]);
## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1,W2,W3), K, T)

## With this unified graph W of size n x n, you can do either spectral clustering or Kernel NMF. If you need help with further clustering, please let us know. 
## for example, spectral clustering

C = 2 					# number of clusters
group = spectralClustering(W, C); 	# the final subtypes information

## you can evaluate the goodness of the obtained clustering results by calculate 
#Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained 
#clustering is very close to the "true" cluster information; 
#if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" 
#cluster information.

displayClusters(W, group);
SNFNMI = calNMI(group, truelabel)

## you can also find the concordance between each individual network and the fused network

ConcordanceMatrix = concordanceNetworkNMI(list(W, W1,W2,W3),2);

################################################################################
# We also provide an example using label propagation to predict 
# the labels of new data points below.
# How to use SNF with multiple views

# Load views into list "dataL"
# load("Digits.RData")
data(Digits)
# Set the other parameters
K = 20 # number of neighbours
alpha = 0.5 # hyperparameter in affinityMatrix
T = 20 # number of iterations of SNF
# Normalize the features in each of the views (optional)
# dataL = lapply(dataL, standardNormalization)

# Calculate the distances for each view
distL = lapply(dataL, function(x) dist2(x, x))

# Construct the similarity graphs
affinityL = lapply(distL, function(x) affinityMatrix(x, K, alpha))
################################################################################
# An example of how to use concordanceNetworkNMI

Concordance_matrix = concordanceNetworkNMI(affinityL, 3);

## The output, Concordance_matrix,  shows the concordance between the fused network and each individual network. 

################################################################################
# Example of how to use SNF to perform subtyping
# Construct the fused network
W = SNF(affinityL, K, T)
# perform clustering on the fused network.
clustering = spectralClustering(W,3);
# use NMI to measure the goodness of the obtained labels.
NMI = calNMI(clustering, label);

################################################################################
# Provide an example of predicting the new labels with label propagation

# Load views into list "dataL" and the cluster assignment into vector "label"
data(Digits)

# Create the training and test data
n = floor(0.8*length(label)) # number of training cases
trainSample = sample.int(length(label), n)
train = lapply(dataL, function(x) x[trainSample, ]) # Use the first 150 samples for training
test = lapply(dataL, function(x) x[-trainSample, ]) # Test the rest of the data set
groups = label[trainSample]

# Set the other
K = 20
alpha = 0.5
t = 20
method = TRUE

# Apply the prediction function to the data
newLabel = groupPredict(train,test,groups,K,alpha,t,method)

# Compare the prediction accuracy
accuracy = sum(label[-trainSample] == newLabel[-c(1:n)])/(length(label) - n)


################################################################################
# References: 
# B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Nature Methods. Online. Jan 26, 2014  
# Website: http://compbio.cs.toronto.edu/SNF/SNF/Software.html

