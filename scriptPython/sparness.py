#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 18:45:14 2017

@author: jsaray
"""
import pandas as pd
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import math
import logging

from syntheticDataFramePreparation import dataFrameCreation
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import Imputer

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

Sc = '[Sparness ] '
S = '[Shapes] - '

# At this point, DataFrames are in variables mRNAToFilter, colToFilter, methColToFilter.

def purgeHighSparsedPatients(df,threshold,barplot=False,title=''):
    logger.debug(Sc+'Purging %s ' %title )
    """Purge patients containing, high sparsed features.
    Patients are Columns.
    Parameters
    ----------
    df : DataFrame
        DataFrame from which we are going to purge high sparsed features
    threshold :  float
        A number between 0 and 1 representing the tolerable sparsity percentage.
        Patients having a percentage of NaN values above the threshold will
        be purged.
    ax : Subplot 
        (optional) a subplot tologger.debug a bar graph representing the 
        sparsity percentage by feature.
        
    Returns
    -------
    numpy array
        A numpy array containing column positions having a NaN percentage higher
        than threshold.
    """
    thr = math.floor(df.shape[0] * threshold)
    colsToDrop = np.array([])
    logger.debug(Sc+'Feature Threshold is %d' % thr)   
    logger.debug(Sc+'Matrix dimensions :  Rows %d , Columns %d'% (df.shape[0],df.shape[1]))
    #axis_x = np.arange(0,df.shape[1])   
    axis_y = np.array([])     
           
    for col in df.columns:
        arr = pd.notnull(df[col])
        nnan = np.sum(arr)      
        axis_y = np.append(axis_y,nnan)
        if (nnan < thr):
            colsToDrop = np.append(colsToDrop,col)
            #df.drop(col,inplace=True,axis=1)
    logger.debug(Sc+'%d patients to drop ' % len(colsToDrop))        
    np.savetxt('debug/sparsePatientsaxis_y.txt',axis_y)
    
    if(colsToDrop == 0):
        logger.debug(Sc+'No patients to drop for %s'% title)
    #ax.bar(axis_x,axis_y)       
    #logger.debug('After purge there are %d columns '% df.shape[1])
    return colsToDrop

def purgeHighSparsedFeatures(df,threshold,barplot=False,title=''):
    """Purge features containing, high sparsed features.
    Features are Rows.
    Parameters
    ----------
    df : DataFrame
        DataFrame from which we are going to purge high sparsed features
    threshold :  float
        A number between 0 and 1 representing the tolerable sparsity percentage.
        Patients having a percentage of NaN values above the threshold will
        be purged.
    ax : Subplot 
        (optional) a subplot tologger.debug a bar graph representing the 
        sparsity percentage by feature.
        
    Returns
    -------
    numpy array
        A numpy array containing column positions having a NaN percentage higher
        than threshold.
    """
    
    thr = math.floor(df.shape[1] * threshold)
    rowsToDrop = np.array([])
    logger.debug(Sc+'Patient Threshold is %d' % thr)   
    logger.debug(Sc+'Matrix dimensions :  Rows %d , Columns %d'% (df.shape[0],df.shape[1]))
    #axis_x = np.arange(0,df.shape[0])   
    axis_y = np.array([])     
    numRows = df.shape[0]       
    for i in range(1,numRows):
        arr = pd.isnull(df.iloc[i])
        nnan = np.sum(arr)  
        axis_y = np.append(axis_y,nnan)
        if (nnan > thr):
            rowsToDrop = np.append(rowsToDrop,i)
    logger.debug ('%d features to drop ' % len(rowsToDrop))
    np.savetxt('debug/sparseFeaturesaxis_y.txt',axis_y)
    #if(barplot):
    #    ax.title.set_text(title)
    #    ax.bar(axis_x,axis_y)       
    #logger.debug('After purge there are %d columns '% df.shape[1])
    return rowsToDrop

def purgeNanEveryWhere(df):
    """Delete rows having NaN everywhere.
    
    Parameters
    ----------
    df : DataFrame
        DataFrame from which rows having NaN everywhere will be deleted.
    
    Returns
    -------
    DataFrame
        The resulting DataFrame after dropping columns and rows containing NaN everywhere.
    """
    #Row-wise dropping
    toDrop = np.array([])
    for i in range(df.shape[0]):
       if( np.sum ( pd.isnull(df.iloc[i]) ) == df.shape[1]-1 ):
           toDrop= np.append(toDrop,i)
    df.drop(df.index[toDrop.astype(int)],inplace=True)       
    #Column-wise dropping
    for col in df.columns:
        arr = pd.notnull(df[col])
        nnan = np.sum(arr)      
        if (nnan == df.shape[1]):
            df.drop(col,inplace=True,axis=1)
    return df

def nanPercentage(df):
    """Compute the ratio of cells having NaN from a DataFrame
    
    Parameters
    ----------
    df: DataFrame
        The original Dataframe
        
    Returns
    -------
    float
        The percentage of cells in NaN
    """
    rows = df.shape[0]
    cols = df.shape[1]
    tot = rows*cols
    
    nanNum = 0
    for i in range(df.shape[0]):
        nanNum = nanNum + np.sum ( pd.isnull(df.iloc[i]) )
    logger.debug ('nan %d tot %d ' % (nanNum, tot) )
    perc = (100*nanNum) / (tot * 1.0)
    return perc

def checkStdDev(df,thr):
    """
    Check that all genes have
    tandard deviation less than a Threshold
    """
    greaterThanThreshold = True
    positions= np.array([])
    for i in range(1,df.shape[0]):
        stdDev = np.std(df.iloc[i,1:].astype(np.longdouble))
        if (stdDev < thr):
            greaterThanThreshold = False
            positions = np.append(positions,i)
    
    return greaterThanThreshold

def computeHist(df,name):
    logger.debug(Sc+'Computing hist for %s' % name)
   
    rows,cols = np.shape(df)
    stds = np.array([])
    logger.debug('getinFor')
    for i in range(1,cols):
        if i % 1000 ==0:
           logger.debug ('Progress '+ str(i) + ' of ' + str(cols))
        arr = df.iloc[1:rows-1,i].astype(float)
        stdev = np.std(arr)
        #logger.debug('stdev %f '% stdev)
        stds = np.append(stds,stdev)
    logger.debug('getOutFor')
    fig = plt.figure()
    plt.hist(stds,bins=40)
    plt.title(name)
    fig.savefig('debug/variableSelection'+name + '.png')

def dropSparsedPatients():
    """
        This function perform the following activities : 
        1. Clean all features And Patients (row and column wise) having NaN everywhere.
        2. RNA table has 0% sparsity. No action is needed.
        3. Purge all patients from miRNA and methylation containing more than a sparsity threshold.
            (This operation implies to erase these patients from the other tables, as we want an intersection)
        4. Purge all features from miRNA and methylation containing more than a sparsity threshold.    
        5. Purge all three dataSets according to a Variance Threshold.
    """
    mRNAToFilter,mirnaToFilter,methColToFilter = dataFrameCreation()


    purgePercentage = 0.6
    logger.debug("Shapes")
    logger.debug(S+ ' %d %d ' % (mRNAToFilter.shape[0],mRNAToFilter.shape[1]))
    logger.debug(S+ ' %d %d ' % (mirnaToFilter.shape[0],mirnaToFilter.shape[1]))
    logger.debug(S+ ' %d %d ' % (methColToFilter.shape[0],methColToFilter.shape[1]))

    logger.info('=========== Purging rows containing only NaN =========\n')
    mirnaToFilter = purgeNanEveryWhere(mirnaToFilter)
    logger.debug("mirna After purging everywhere %d %d " % (mirnaToFilter.shape[0],mirnaToFilter.shape[1]))
    methColToFilter = purgeNanEveryWhere(methColToFilter)
    logger.debug("meth After purging everywhere %d %d " % (methColToFilter.shape[0],methColToFilter.shape[1]))
    
    logger.debug('miRNA After purging NaN column and rowise %d %d '% (mirnaToFilter.shape[0],mirnaToFilter.shape[1]))
    logger.debug('Methylation After purging NaN column and rowise %d %d '% (methColToFilter.shape[0],methColToFilter.shape[1]))  
    
    logger.info('\n=========== Purging Patients in RNA table above %f of Sparsity\n' % purgePercentage )
    colsMrna = purgeHighSparsedPatients(mRNAToFilter,purgePercentage,title='rna')
    if len(colsMrna) != 0:
        mRNAToFilter.drop(colsMrna,inplace=True,axis=1)
    # We know already that the sparsity percentage for this Matrix is 0 so no need to filter out
    # any patient in other tables, however a general approach should contain code
    # filtering out patients on other tables.
    
    logger.debug('=========== Purging Patients in miRNA table above %f of Sparsity \n' % purgePercentage ) 
    
    cols = purgeHighSparsedPatients(mirnaToFilter,purgePercentage,title='mirna')
    rows = purgeHighSparsedFeatures(mirnaToFilter,purgePercentage,title='mirna')

    logger.debug('High Sparsed Patients %d' % len(cols))
    logger.debug('High Sparsed Features %d' % len(rows) )   
    mi_to_rnaCols = np.array([])
    mi_to_methCols = np.array([])
        
    for col in cols:
        mi_to_rnaCols = np.append(mi_to_rnaCols,col[:12])
        mi_to_methCols = np.append(mi_to_methCols,col[:12])
    
    # Drop miRNA sparse patients in RNA table        
    mRNAToFilter.drop(mi_to_rnaCols,inplace=True,axis=1)
    logger.debug('final dim rna ')
    logger.debug(mRNAToFilter.shape)
    mirnaToFilter.drop(cols,inplace=True,axis=1)
    mirnaToFilter.drop(mirnaToFilter.index[rows.astype(int)],inplace=True) 
    logger.info('final dim mirna')
    logger.debug(mirnaToFilter.shape)
    logger.debug('methcols')
    np.savetxt('debug/methcols_str.txt',methColToFilter.columns,delimiter=" ", fmt="%s")
    logger.debug(methColToFilter.columns)
    # Drop miRNA sparse patients in Methylation table.
    methColToFilter.drop(mi_to_methCols,inplace=True,axis=1)
    methColToFilter.to_csv('debug/afterMEth.csv')
    logger.debug('final dim meth')
    logger.debug(methColToFilter.shape)
    logger.debug('=========== Purging Patients in Methylation table above %f of Sparsity ====== \n' % purgePercentage )  
    colsMeth = purgeHighSparsedPatients(methColToFilter,purgePercentage,title='meth')
       
    rowsMeth = purgeHighSparsedFeatures(methColToFilter,purgePercentage,title='mirna')
    
    logger.debug('High Sparsed Patients %d' % len(colsMeth))
    logger.debug('High Sparsed Features %d' % len(rowsMeth) )
    
    if(len(colsMeth)!= 0):
        methColToFilter.drop(colsMeth,inplace=True,axis=1)
    if(len(rowsMeth)!= 0):
        methColToFilter.drop(methColToFilter.index[rowsMeth.astype(int)],inplace=True) 
    
    
    logger.debug('\n After cleaning Dimensions \n=========================')
    logger.debug('After Cleaning rna  %d %d' % (mRNAToFilter.shape[0],mRNAToFilter.shape[1]))
    logger.debug('After Cleaning mirna %d %d' % (mirnaToFilter.shape[0],mirnaToFilter.shape[1]))
    logger.debug('After Cleaning meth %d %d' % (methColToFilter.shape[0],methColToFilter.shape[1]))
    logger.debug('==========================\n')
     # We know already that for a 0.6 sparsity percentage we dont need to filter out
    # any patient, so no need to do and drop patients on other tables, however a general approach should contain code
    # filtering out patients on other tables.
    
    
    mRNAToFilter = mRNAToFilter.T
    mRNAToFilter.columns = mRNAToFilter.iloc[0,:].axes
    mRNAToFilter.set_index( mRNAToFilter.iloc[:,0]) 
    
    logger.debug('Transposed Shape %d %d' % (mRNAToFilter.shape[0],mRNAToFilter.shape[1]))
    logger.debug('First column ')
    logger.debug(mRNAToFilter.iloc[:,0].axes)
    
    rowsRna,colsRna = mRNAToFilter.shape
    logger.debug('just before hist %d %d ' % (rowsRna,colsRna))
    computeHist(mRNAToFilter.iloc[1:rowsRna, 0:colsRna].astype(float),'StdDev RNA')
    mRNAToFilter.to_csv('debug/debugrna.csv')
    
    logger.info('Variance Threshold filtering ( methylation is already preprocessed in Firebrowse)')
    vt = VarianceThreshold(threshold=0.5)
    
    
    logger.debug(S+' Number of cols in RNA before transfo %d' % mRNAToFilter.shape[1])
    
   
    filteredRna = vt.fit_transform(mRNAToFilter.iloc[1:rowsRna, 0:colsRna].astype(float))
    labelsRna = [mRNAToFilter.columns[i] for i in np.arange(0,len(vt.get_support())) if vt.get_support()[i]]  
     
    logger.debug('The Variances ARE ')
    logger.debug(S+' Number of cols after RNA transfo %d' % filteredRna.shape[1])
    ##### OJOOOOOOOOOOOO NO SE DEBERIA HACER IMPUTING AQUI ##########################
    ##### TODOOOOOOOOOOOOOOOOO ############################
    imp = Imputer()
    mirnaToFilter = mirnaToFilter.T
    mirnaToFilter.columns = mirnaToFilter.iloc[0,:]
    mirnaToFilter.set_index( mirnaToFilter.iloc[:,0]) 
    
    rowsMiRna,colsMiRna = mirnaToFilter.shape
    imputedColToFilter = imp.fit_transform(mirnaToFilter.iloc[1:rowsMiRna , 0:colsMiRna].astype(float))
    computeHist(pd.DataFrame(data=imputedColToFilter[1:rowsMiRna , 1:colsMiRna].astype(float)), 'StdDev miRNA')

    logger.debug(S+' Number of cols in miRNA before transfo %d' % colsMiRna)
    
    
   
    
    filteredMiRna = vt.fit_transform(imputedColToFilter)
    labelsmiRNA = [mirnaToFilter.columns[i] for i in np.arange(0,len(vt.get_support())) if vt.get_support()[i]] 
    
    
    logger.debug(S+' Number of cols in miRNA after transfo %d' % filteredMiRna.shape[1])
    
    methColToFilter = methColToFilter.T
    methColToFilter.columns = methColToFilter.iloc[0,:]
    methColToFilter.set_index( methColToFilter.iloc[:,0]) 
    
    rowsMeth,colsMeth = methColToFilter.shape
    logger.debug(' %d %d just before transform meth '%(rowsMeth,colsMeth))
    imputedMethColToFilter = imp.fit_transform(methColToFilter.iloc[1:rowsMeth, 0:colsMeth].astype(float))

    logger.debug('--dim rna %d %d' % (filteredRna.shape[0],filteredRna.shape[1])   )
    
    logger.debug('--dim mirna %d %d'% (filteredMiRna.shape[0],filteredMiRna.shape[1])  )

    logger.debug('--dim meth %d %d'% (imputedMethColToFilter.shape[0],imputedMethColToFilter.shape[1]))
    
    
    logger.debug(" list comprehension ")

   

    METH_COLS  = methColToFilter.iloc[0,:].values.copy()
    
    logger.debug("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO TO WRITE COLUMNS 0000000000000000000000000000000")
    logger.debug(len(labelsRna))     
    logger.debug(len(METH_COLS))   
    logger.debug(len(labelsmiRNA))   

    #percRna = nanPercentage(mRNAToFilter)
    ##percMiRna = nanPercentage(colToFilter)
    #percMeth = nanPercentage(methColToFilter)
            
    #print ('Sparsity Percentages : ')
    #print ('mRNA %f ' % percRna )
    #print ('miRNA %f ' % percMiRna )
    #print ('percMeth %f ' % percMeth )
    #print ( " ======== Check all Methylation probes have StdDev > 0.2)======== " )
    #std = checkStdDev(methColToFilter,0.2)
        
    #if(std):
    #   logger.debug ('Threshold true for all genes')
    #else:
    #   logger.debug ('Threshols false for all genes')
    rnaDf_toWrite = pd.DataFrame(data=filteredRna,columns=labelsRna,index=mRNAToFilter.index.values[1:])
    mirnaDf_toWrite = pd.DataFrame(data=filteredMiRna,columns=labelsmiRNA,index=mirnaToFilter.index.values[1:])
    methDf_toWrite = pd.DataFrame(data=imputedMethColToFilter,columns=METH_COLS,index=methColToFilter.index.values[1:])
    
    
    rnaDf_toWrite=rnaDf_toWrite.sort_index()  
    mirnaDf_toWrite=mirnaDf_toWrite.sort_index()
    methDf_toWrite=methDf_toWrite.sort_index()     
    #fig.set_size_inches(20, 20)
    #fig.tight_layout()
    #fig.savefig('myfig.png')
    
    logger.debug('mrna Columns')
    
    rna_toWrite = pd.DataFrame(data=rnaDf_toWrite)
    mirna_toWrite = pd.DataFrame(data=mirnaDf_toWrite)
    meth_toWrite = pd.DataFrame(data=methDf_toWrite)
    
    rna_toWrite.to_csv('debug/readyToMixRNA.csv')
    mirna_toWrite.to_csv('debug/readyToMixMIRNA.csv')
    meth_toWrite.to_csv('debug/readyToMixMETH.csv')

    return (mRNAToFilter,mirnaToFilter,methColToFilter)
dropSparsedPatients()