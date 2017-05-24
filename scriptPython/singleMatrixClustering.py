#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:59:53 2017

@author: jsaray
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging

from mpl_toolkits.mplot3d import Axes3D

from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans


from sklearn.metrics.pairwise import euclidean_distances


from sklearn import metrics
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from pyclustering.cluster import xmeans

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def gaussianMixtures(nClusters):
    iris = datasets.load_iris()
    X = iris.data
    gmm = GaussianMixture(n_components=3, covariance_type='spherical').fit(X)
    labels = gmm.predict(X)
    probas = gmm.predict_proba(X)
    fig = plt.figure(1, figsize=(8, 6))
    X_reduced = PCA(n_components=3).fit_transform(X)
    nlabels = list(map(change,labels))
    logger.debug(nlabels)
    ax = Axes3D(fig, elev=-150, azim=110)
    ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=nlabels,cmap=plt.cm.Paired)
    ax.set_title("First three PCA directions")
    ax.set_xlabel("1st eigenvector")
    ax.w_xaxis.set_ticklabels([])
    ax.set_ylabel("2nd eigenvector")
    ax.w_yaxis.set_ticklabels([])
    ax.set_zlabel("3rd eigenvector")
    ax.w_zaxis.set_ticklabels([])
    plt.show()
    logger.debug(type(X))
    logger.debug(X)
    logger.debug (labels)
    logger.debug ('Probas are')
    logger.debug (probas)
    np.savetxt('debug/probas.txt',probas)

def computeJaccardMatrix(dict1,dict2):
    cardunion = np.union1d(dict1,dict2).size()
    cardinter = np.intersect1d(dict1,dict2).size()
    jaccard = cardinter*1.0 / cardunion
    return jaccard

def change(x):
    if (x == 0):
        ret = 'm'
    if (x == 1):
        ret = 'y'
    if (x == 2):
        ret = 'b'
    if (x == 3):
        ret = 'r'
    return ret

def findOptimalClusterNumber(df,method,algo,omics):
    logger.debug('Shape of DF %d %d' % (np.shape(df)[0],np.shape(df)[1]) )
    sils = np.array([])
    elbow = np.array([])
    if(method == 'silhouette'):
        if (algo == 'spectral'):
            for k in np.arange( 2, 100) :
                logger.debug('Computing for %d Clusters ' % k)
                spect = SpectralClustering(n_clusters=int(k), eigen_solver='arpack', affinity="nearest_neighbors")
                clusters = spect.fit(df)
                silhouette = metrics.silhouette_score(df,clusters.labels_)
                sils = np.append(sils,silhouette)
                nam = 'debug/silhouette_'+algo+'-'+method+'-'+omics + '.txt'
                np.savetxt(nam,sils)
        if (algo == 'agglom'):
            for k in np.arange( 2, 100) :
                logger.debug('Computing for %d Clusters ' % k)
                spect = hierarchicalClustering(df,'Agglomerative',k)
                
                silhouette = metrics.silhouette_score(df,spect)
                sils = np.append(sils,silhouette)
                nam = 'debug/silhouette_'+algo+'-'+method+'-'+omics +'.txt'
                np.savetxt(nam,sils)
                
    if(method == 'gap'):
        if (algo == 'spectral'):
            for k in np.arange( 2,100) :
                logger.debug('Computing for %d Clusters ' % k)
                spect = SpectralClustering(n_clusters=k, eigen_solver='arpack', affinity="nearest_neighbors")
                clusters = spect.fit(df)
                silhouette = metrics.silhouette_score(df,clusters.labels_)
                sils = np.append(sils,silhouette)
                np.savetxt('debug/silhouette.txt',sils)
    if(method == 'elbow'):
        if (algo == 'kmeans'):
            for k in np.arange( 2,100) :
                logger.debug('Computing for %d Clusters ' % k)
                spect = KMeans(n_clusters=int(k),n_init=3)
                clusters = spect.fit(df)
                el = clusters.inertia_
                elbow = np.append(elbow,el)
                np.savetxt('debug/elbow_' + method + '_' +algo + '_' +omics + '.txt',elbow)
    

def kMeansClustering(df,title,nClusters):
    """Performs kmeans clustering
    
    Parameters
    ----------
    df : DataFrame
        DataFrame we use in our clustering algorithm
    title : string
        Our figure title
    """
    plt.title(title)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    spect = KMeans(n_clusters=nClusters, n_init=10)
    logger.debug('Invoking Spectral clustering on a matrix having %d, %d ' % (df.shape[0],df.shape[1]))
    clusters = spect.fit(df)
    logger.debug('Result SpectralClustering')
    logger.debug(clusters.labels_)
    xs=df[:,0]
    ys=df[:,1]
    zs=df[:,2]
    ax.scatter(xs, ys, zs, c=clusters.labels_,marker='o')
    fig.savefig('debug/' + title + '.png')
    return clusters.labels_

def spectralClustering(df,title,nClusters):
    """Perform Clustering using Spectral Clustering technique.
    
    Parameters
    ----------
    df : DataFrame
        DataFrame we use in our clustering algorithm.
    title : string 
        Our figure title.      
    """
    
    spect = SpectralClustering(n_clusters=nClusters, eigen_solver='arpack', affinity="nearest_neighbors")
    logger.debug('Invoking Spectral clustering on a matrix having %d, %d ' % (df.shape[0],df.shape[1]))
    clusters = spect.fit(df)
    logger.debug('Result SpectralClustering')
    logger.debug(clusters.labels_)

    return clusters.labels_

def hierarchicalClustering(df,title,nClusters):
    """Perform Clustering using hierarchical Clustering technique.
    
    Parameters
    ----------
    df : DataFrame
        DataFrame we use in our clustering algorithm.
    title : string 
        Our figure title.      
    """
    #plt.title(title)
    #fig = plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    aggloClustering = AgglomerativeClustering(n_clusters=nClusters,affinity='euclidean',linkage='average')
    clusters = aggloClustering.fit(df)
    logger.debug('Invoking HIerarchical - Agglo clustering on a matrix having %d, %d ' % (df.shape[0],df.shape[1]))
    logger.debug('Result HIerarchical - Agllo Clustering')
    logger.debug(clusters.labels_)
    #xs=df[:,0]
    #ys=df[:,1]
    #zs=df[:,2]
    #ax.scatter(xs, ys, zs, c=clusters.labels_,marker='o')
    #fig.savefig('debug/' + title + '.png')
    return clusters.labels_
    
def densityClustering(df,title):
    """Perform Clustering using density Clustering technique.
    
    Parameters
    ----------
    df : DataFrame
        DataFrame we use in our clustering algorithm.
    title : string 
        Our figure title.  
    """
    
    densityClustering = DBSCAN(120, min_samples=20)
    clusters = densityClustering.fit_predict(df)
    logger.debug('Invoking Spectral clustering on a matrix having %d, %d ' % (df.shape[0],df.shape[1]))

   
    return clusters

THRESHOLD = 0.6
CLUSTER_NUMBER = 6
mRNA = pd.read_csv('/Users/jsaray/DataFrameCreation/readyToMixRNA.csv',sep=',')
(rnaRows,rnaCols) = mRNA.shape
mRNA = mRNA.iloc[0:rnaRows,1:rnaCols]

MIRna = pd.read_csv('/Users/jsaray/DataFrameCreation/readyToMixMIRNA.csv',sep=',')
(mirnaRows,mirnaCols) = MIRna.shape
MIRna = MIRna.iloc[0:mirnaRows,1:mirnaCols]

meth = pd.read_csv('/Users/jsaray/DataFrameCreation/readyToMixMETH.csv',sep=',')
(methRows,methCols) = meth.shape
meth = meth.iloc[0:methRows,1:methCols]

logger.debug('BEFORE FILTERING ')

logger.debug('Dimensions rna ' )
logger.debug( mRNA.shape)
logger.debug('Dimensions mirna ' )
logger.debug( MIRna.shape)
logger.debug('Dimensions meth ' )
logger.debug( meth.shape)

#subMatrRna = mRNA.iloc[0:mRNA.shape[0],1:mRNA.shape[1]]

#matr = euclidean_distances(subMatrRna)
#np.savetxt('debug/eumatr.csv',matr,delimiter=',')
#clus = densityClustering(subMatrRna,'rna')
#logger.debug(len(set(clus)))

#clus = densityClustering(subMatrRna,'rna')
#logger.debug(len(set(clus)))

#clus = densityClustering(subMatrRna,'rna')
#ogger.debug(len(set(clus)))


#findOptimalClusterNumber(mRNA,'silhouette','spectral','rna')
#findOptimalClusterNumber(mRNA,'silhouette','agglom','rna')

#findOptimalClusterNumber(MIRna,'silhouette','spectral','mirna')
#findOptimalClusterNumber(MIRna,'silhouette','agglom','mirna')

#findOptimalClusterNumber(meth,'silhouette','spectral','meth')
#findOptimalClusterNumber(meth,'silhouette','agglom','meth')


findOptimalClusterNumber(mRNA,'elbow','kmeans','rna')
findOptimalClusterNumber(MIRna,'elbow','kmeans','mirna')
findOptimalClusterNumber(meth,'elbow','kmeans','meth')
'''
logger.debug(' Spectral Clustering ')
sp_RNA_labels = spectralClustering(filteredRna,'Spectral RNA',CLUSTER_NUMBER)
logger.debug('rna label valuees')
logger.debug(sp_RNA_labels)
sp_miRNA_labels = spectralClustering(filteredMiRna,'Spectral miRNA',CLUSTER_NUMBER)
sp_meth_labels = spectralClustering(imputedMeth,'Spectral Meth',CLUSTER_NUMBER)

for i in np.arange(0,CLUSTER_NUMBER):
    mask = np.where(sp_RNA_labels == i)

jaccMatrix = computeJaccardMatrix(sp_RNA_labels)



hier_RNA_labels = hierarchicalClustering(filteredRna,'Hierarchical - RNA',CLUSTER_NUMBER)
hier_MIRNA_labels = hierarchicalClustering(filteredMiRna,'Hierarchical - miRNA',CLUSTER_NUMBER)
hier_METH_labels = hierarchicalClustering(imputedMeth,'HIerarchical - Meth',CLUSTER_NUMBER)

dens_RNA_labels = densityClustering(filteredRna,'Density - RNA')
dens_miRNA_labels = densityClustering(filteredMiRna,'Density - miRNA')
dens_Meth_labels = densityClustering(imputedMeth,'Density - Meth')'''


