#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 16:22:53 2017

@author: jsaray
"""

from dataFrameCreation import dataFrameCreation
from sklearn.decomposition import PCA
import pandas as pd 

def fusion(rna,mirna,meth):
    
    rna = pd.DataFrame.transpose(rna)
    mirna = pd.DataFrame.transpose(mirna)
    meth = pd.DataFrame.transpose(meth)
    
    frames = [rna,mirna,meth]
    fusionM = pd.concat(frames,axis=1)
    return fusionM


dfRna, dfMiRna, dfMeth = dataFrameCreation()
print ('Features rna, mirna, meth %d %d %d ' % (dfRna.shape[1],dfMiRna.shape[1],dfMeth.shape[1]) )
bigMatrix = fusion(dfRna,dfMiRna,dfMeth)
print('Fusion Dimensions ')
print(bigMatrix.shape)