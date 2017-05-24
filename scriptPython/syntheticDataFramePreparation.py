#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 12:11:57 2017

@author: jsaray
"""
import pandas as pd
import numpy as np
import logging
from intersectionComputing import getIntersection
from collections import Counter

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

Sc = '[SyntheticDataFramePreparation] '
S = '[Shapes] - '


def checkForMultiSamplePatients(patients):
    """
    This function retrieves a dictionary with all the patients having more
    than one sample in the DataSet
    Parameters
    ----------
    patients : array
        Patients array to compute occurrences for each patient
    
    Returns
    -------
    dictionary
        Dictionary containing all patients occuring more than once in the array.
    """
    substrArray = [c[:12] for c in patients]
    # Counter is a nice Python class that counts number of occurrences in an array
    count = Counter(substrArray)

    # New C is the number of patients having more than one occurence in DataSet,
    # this value NEVER must be different than zero, if it is, it means there are
    # some patients having more than one sample, we have then to choose only one sample 
    # per patient.
    newc = {k: v for k, v in count.items() if v > 1}
    return newc

def leaveOneColumnPerPatient(df,dic):
    """
    Some data frames contain patients having more than one Sample.
    This function guarantees we have ONE and ONLY ONE sample per patient.
    Parameters
    ----------
    df : DataFrame
        data frame to analyze
    dic : Dictionary
        Dictionary containing which patients are occurring more than once.
    
    Returns
    -------
    DataFrame
        A dataframe containing only one sample per patient.
    """
    logger.debug('just get in %s'%dic)
    logger.debug('Data Frame dimensions %s %s ' % (df.shape[0],df.shape[1]))
    for k,v in dic.items():
        recreat = [c for c in df.columns if k in c]
        diffPos = np.array([])
        for c in recreat:
            pos = df.columns.get_loc(c)
            diffPos = np.append(diffPos,pos)
        df.drop(df.columns[diffPos[1:].astype(int)],inplace=True,axis=1)
    logger.debug('final length after droppin %d' % len(df.columns))    
    return df

def dataFrameCreation():
    """
    This function creates three DataFrames (rna, miRna and Methylation)
    
    Returns
    -------
    Tuple
        A tuple of three data frames.
    """
    intersection = getIntersection()
    intersectionRna = np.insert(intersection,0,'feature')
    logger.debug('DEEEEEEBUGGGGGG')
    logger.info('Intersection done. ')
    logger.info('Intersection length. %d'% len(intersection))
    
    # mRNA matrix N x Mk1
    logger.info('======================== BUILDING RNA MATRIX ====================')
    mRNA = pd.read_csv('/Users/jsaray/DataFrameCreation/TCGACRC_expression-merged.tsv',sep='\t')
    mRNA.index = mRNA['feature']
    logger.debug('Before Dropping %d %d '%(mRNA.shape[0],mRNA.shape[1]))
    mRNAToFilter = mRNA.filter(intersectionRna,axis=1)
    logger.debug(Sc + 'After Dropping %d %d '% (mRNAToFilter.shape[0],mRNAToFilter.shape[1]))
    logger.debug(Sc +S+'In memory - RNA Dimensions %d %d' % (mRNAToFilter.shape[0],mRNAToFilter.shape[1]))
    #miRNA matrix N x Mk2
    logger.info(Sc + '======================== BUILDING MIRNA MATRIX , (log transformed)=== ')
    MIRna = pd.read_csv('/Users/jsaray/DataFrameCreation/miRseq_mature_RPM_log2.csv',sep='\t')
    MIRna.index = MIRna['feature']
    
    
    logger.debug(S+'Before Dropping %d %d  '%(MIRna.shape[0],MIRna.shape[1]))
    
    # fixed is the size of patients contained in the intersection AND in miRNA
    intersect = [c for c in MIRna.columns if c[:12] in intersection]
    newc = checkForMultiSamplePatients(intersect)
    logger.debug('NUMBER OF PATIENTS REPEATED IN MIRNA : ')
    logger.debug(len(newc))
    
    logger.debug('Columns of miRNA before filtering')
    
    colToFilter = MIRna.ix[:, intersect]
    logger.debug('col to filter columns ')
    logger.debug(str(len(colToFilter.columns)))
    
    logger.debug('calling leave one column')
    colToFilter = leaveOneColumnPerPatient(colToFilter,newc)   
    colToFilter.insert(0,'feature',MIRna.iloc[:,0])
    
    logger.debug(S+'After Dropping %d %d'% (colToFilter.shape[0],colToFilter.shape[1]))
    
    logger.info('====================== BUILDING METHYLATION MATRIX ==================')
    meth = pd.read_csv('/Users/jsaray/DataFrameCreation/COADREAD.meth.by_max_stddev.data.csv',sep='\t')
    meth.index = meth['feature']
    logger.debug(Sc + ' before dropping columns %d %d ' %(meth.shape[0],meth.shape[1]))
    logger.debug(' Renaming Methylation columns ' )
   
    for i in range(0,meth.shape[0]):
        meth.iloc[i,0] = 'Beta'+ meth.iloc[i,0] 
   
    
    logger.debug(Sc + 'before filtering')
    logger.debug(meth.shape)
    # fixed is the size of patients contained in the intersection AND in miRNA
    fixed = [c for c in meth.columns if c[:12] in intersection]
    np.savetxt('fixed.txt',fixed,delimiter=" ",fmt="%s")
    newc = checkForMultiSamplePatients(fixed)
    
    logger.debug(Sc + 'NUMBER OF DUPLICATES FOR METHYLATION')
    logger.debug(Sc + '%d' % len(newc))
    
    methColToFilter = meth.ix[:,fixed]
    logger.debug('After filtering dimensions %d %d' % (methColToFilter.shape[0],methColToFilter.shape[1]))
    np.savetxt('justBefore.txt',methColToFilter.columns,delimiter=" ",fmt="%s")
    methColToFilter = leaveOneColumnPerPatient(methColToFilter,newc)
    np.savetxt('afterLeaveOneColumn.txt',methColToFilter.columns,delimiter=" ",fmt="%s")
    methColToFilter.insert(0,'feature',meth.iloc[:,0])
    methColToFilter.to_csv('debug/methFiltered.csv')
    logger.debug('Filtered data frame dimensions : ')
    logger.debug(Sc)
    logger.debug(methColToFilter.shape)
    
    return mRNAToFilter,colToFilter,methColToFilter