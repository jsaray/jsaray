#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:01:29 2017

@author: jsaray
"""

from sparness import dropSparsedPatients
from sklearn.decomposition import PCA
from sklearn.preprocessing import Imputer
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

rna,mirna,meth = dropSparsedPatients()

rna.to_csv('debug/rnaFirstCheck.csv')
mirna.to_csv('debug/mirnaFirstCheck.csv')
meth.to_csv('debug/methFirstCheck.csv')



print('shape rna')
print(rna.shape)
print('shape mirna')
print(mirna.shape)
print('shape meth')
print(meth.shape)


print('rna dim')
print(rna.shape)

rna = rna.T
print('transp rna dim')
print(rna.shape)
print('rna columnsssss')
print(rna.columns)
print('rna values keys')
print(rna.index.values)

rna.to_csv('debug/rnaT.csv')

rnaToMix = rna
rnaToMix.to_csv('debug/rnaAfterDrop.csv')

mirna = mirna.T

mirnaToMix = mirna

mirnaToMix.to_csv('debug/mirnaAfterDrop.csv')
#mirna.columns = mirna.iloc[0,:]

meth = meth.T
methToMix = meth

methToMix.to_csv('debug/methAfterDrop.csv')

print('rna cols')
print(rna.columns)
print('mirna cols')
print(mirna.columns)
print('meth cols')
print(meth.columns)

frames = [rnaToMix,mirnaToMix,methToMix]
bigdf = pd.concat(frames,axis=1)



print('bigdf dimensions')
print(bigdf.shape)
pca = PCA(n_components=0.7)

bigdf.to_csv('debug/CAMERA_READY_HEADER_LAST_ROW.csv')

lastline = bigdf.iloc[bigdf.shape[0] - 1]
newIndexes = np.arange(1,bigdf.shape[0] )
newIndexes = np.append(newIndexes,0)
bigdf.index = newIndexes
bigdf.drop(bigdf.index[bigdf.shape[0] - 1])

bigdf.to_csv('debug/CAMERA_READY.csv')

numlin = bigdf.shape[0] - 1
pcaMatrix = bigdf.iloc[0:numlin,1:]
matr = pcaMatrix.values


stdevs = np.std(matr)

print('type of mtr')
print(type(matr))
print('dimension')
print(np.shape(matr))

imp = Imputer()

impPcaMatrix = imp.fit_transform(matr)

standardScaler = StandardScaler()

impPcaMatrix = standardScaler.fit_transform(impPcaMatrix)

print('imputer dim')
print(impPcaMatrix.shape)
transfo  = pca.fit_transform(impPcaMatrix)
print('Transformation type ')
print(type(transfo))
print('dimensions')
print(np.shape(transfo))
print('content')
print(transfo)


print('Variance ratio')
print(pca.explained_variance_ratio_)
print('Number of components')
print(pca.n_components)

xs = transfo[:,0]
ys = transfo[:,1]
zs = transfo[:,2]

subm = np.vstack((xs,ys,zs))
print('dimsubm')
print(np.shape(subm))
subm = np.transpose(subm)

Kmeans = KMeans(n_clusters = 4,random_state=0).fit(subm)

print(Kmeans.labels_)

ax.scatter(xs, ys, zs, c=Kmeans.labels_,marker='o')

fig.savefig('debug/clustering.png')