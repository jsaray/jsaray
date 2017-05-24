#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 17:13:57 2017

@author: jsaray
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import Imputer

from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN


def spectralClustering(df,title):
    plt.title(title)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    spect = SpectralClustering(n_clusters=4, eigen_solver='arpack', affinity="nearest_neighbors")
    clusters = spect.fit(df)
    print('Result SpectralClustering')
    print(clusters.labels_)
    plt.plot(clud)
    ax.scatter(xs, ys, zs, c=clusters.labels_,marker='o')
    fig.savefig('debug/clustering.png')
    
def hierarchicalClustering(df):
    hierch = AgglomerativeClustering()
    
def dbScan():
    

def stdInThread(matr,i):
    return np.std(matr[:,i])

def computeHist(df,name):
    print('Computing hist for %s' % name)
   
    rows,cols = np.shape(df)
    stds = np.array([])
    print('getinFor')
    for i in range(1,cols):
        if i % 1000 ==0:
            print ('Progress '+ str(i) + ' of ' + str(cols))
        arr = df.ix[1:rows-1,i].astype(float)
        stdev = np.std(arr)
        #print('stdev %f '% stdev)
        stds = np.append(stds,stdev)
    print('getOutFor')
    fig = plt.figure()
    plt.hist(stds,bins=40)
    plt.title(name)
    fig.savefig('debug/variableSelection'+name + '.png')
    
   

rnaToMix = pd.read_csv('debug/rnaAfterDrop.csv')
print('dim rna ')
print(np.shape(rnaToMix))
mirnaToMix = pd.read_csv('debug/mirnaAfterDrop.csv')
print('dim mirna ')
print(np.shape(mirnaToMix))
methToMix  = pd.read_csv('debug/methAfterDrop.csv')
print('dim meth ')
print(np.shape(methToMix))

computeHist(methToMix,'METH') 

vtRna = VarianceThreshold(threshold=0.5)
vtMiRna = VarianceThreshold(threshold=0.5)
vtMeth = VarianceThreshold(threshold=0.05)

rowsRna = rnaToMix.shape[0]
colsRna = rnaToMix.shape[1]

filteredRna = vtRna.fit_transform(rnaToMix.ix[1:rowsRna-1, 1:colsRna].astype(float))

rowsMiRna = mirnaToMix.shape[0]
colsMiRna= mirnaToMix.shape[1]

imp = Imputer()
imputedMiRna = imp.fit_transform(mirnaToMix.ix[1:rowsMiRna-1 , 1:colsMiRna].astype(float))
filteredMiRna = vtMiRna.fit_transform(imputedMiRna)

rowsMeth = methToMix.shape[0]
colsMeth = methToMix.shape[1]

imp2 = Imputer()
imputedMeth= imp2.fit_transform( methToMix.ix[1:rowsMeth-1, 1:colsMeth].astype(float) )

filteredMeth = vtMeth.fit_transform( imputedMeth )

frames = [filteredRna,filteredMiRna,filteredMeth]
allTogether = pd.concat(frames,axis=1)

print ('after rna')
print(np.shape(filteredRna))
print ('after mirna')
print(np.shape(filteredMiRna))
print ('after meth')
print(np.shape(filteredMeth))

spectralClustering(filteredRna,'Spectral Clustering, RNA')
spectralClustering(filteredMiRna,'Spectral Clustering, miRNA')
spectralClustering(filteredMeth,'Spectral Clustering, Methylation')
spectralClustering(allTogether,'Spectral Clustering, all Together')


           
           
           
