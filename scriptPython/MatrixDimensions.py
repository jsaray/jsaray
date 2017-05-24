#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 12:30:29 2017

@author: jsaray
"""
import pandas as pd
import numpy as np

def number(cad):
    arr = cad
    new = np.array([])
    print ('tot ' +str(len(arr)))
    print(arr)
    for s in arr:
        new = np.append(new, s[8:12])
    
    print ('new is')
    print (new)
    un = np.unique(new)
    print('different values ' + str(len(un)))
    

df = pd.read_csv('/Users/jsaray/Internship/Code/data/dataSamples/RNASeq/GenesData/Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/COADREAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.csv',sep='\t',header=None)

arr = df.iloc[0].as_matrix()
print(arr)

print (type(arr))
print (arr.size)
number(arr[1:])