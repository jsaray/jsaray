#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 16:35:07 2017

@author: jsaray
"""
import numpy as np
import re
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def firstLineOfFile(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    return first_line

def number(cad):
    """Given a string containing a series of TCGA barcodes, return
    an array with each TCGA barcode appearing only once.
    
    Parameters
    ----------
    cad : string containing different tcga barcodes.
    
    Returns
    -------
    Array
        An array containing each TCGA barcode only once.
    """
    arr = re.findall(r"[\w-]+",cad)
    new = np.array([])
    logger.debug('Tot ' +str(len(arr)))

    for s in arr:
        new = np.append(new, s[:12])
    
    logger.debug('Four digits array')
    un = np.unique(new)
    logger.debug('Different values ' + str(len(un)))
    return un

def getIntersection():
    """
    Computes the patient intersection between 3 DataSets. That means, patients
    that are present on 3 DataSets.
    
    Returns
    -------
    Array
        An array containing intersection of 3 DataSets.
    """
    logging.info('========= miRNA merging =============')
    FILE_PREFIX = '/Users/jsaray/DataFrameCreation/'
    miRNA = firstLineOfFile(FILE_PREFIX + 'miRseq_mature_RPM_log2.csv' )

    miRnaMerged = number(miRNA)

    logging.info('=========== mRna already merged ==========')

    mRNAMerged = firstLineOfFile(FILE_PREFIX + 'TCGACRC_expression-merged.tsv')
    mRNAMerged = number(mRNAMerged)
    intermiRNAmRNA = np.intersect1d(mRNAMerged,miRnaMerged)
    logging.info('len of inter miRNA mRNA ' + str(len(intermiRNAmRNA)))
    
    logging.info('================== Methylation ============')
    methMerged = firstLineOfFile(FILE_PREFIX + 'COADREAD.meth.by_max_stddev.data.csv')
    methMerged = number(methMerged)
    allMerged = np.intersect1d(intermiRNAmRNA,methMerged)
    logging.info('unique length')
    logging.info(len(np.unique(allMerged)))
    allMerged =  allMerged[1:len(allMerged)-1]
    logging.info('len of all Merged ' + str(len(allMerged)))
    np.savetxt('debug/inter.txt',allMerged,delimiter=" ",fmt="%s")
    np.savetxt('debug/meth.txt',methMerged,delimiter=" ",fmt="%s")
    return allMerged
getIntersection()