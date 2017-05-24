#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:10:44 2017

@author: jsaray
"""
#https://gdc-api.nci.nih.gov/cases?facets=files.experimental_strategy&from=1&size=0&pretty=true
import pandas as pd
import firebrowse
import requests
import string
import json
from firebrowse import fbget
import numpy as np
import urllib2

def createFile(filename,strin):
    output = open('/Users/jsaray/Internship/Code/'+filename,'w')
    output.write(strin)
    output.close()
    
def createFileFromArray(filename,arr):
    output = open('/Users/jsaray/Internship/Code/'+filename,'w')
    for line in arr:
        output.write(line+'\n')
    output.close()

def filterAnswer(pattern,Text):
    toRet = []
    stri = string.split(Text,'\n')

    print (str(len(stri)))
    for line in stri:
        if pattern in line.lower(): 
            toRet.append(line)
    return toRet
cases_endpt = 'https://gdc-api.nci.nih.gov/files'
BIGPRX = "http://www.cbioportal.org/webservice.do?cmd="
'''
print ('Methylation FIREBROWSE')
A = firebrowse.Archives().StandardData(format='json',data_type='Methylation',cohort='COADREAD')
createFile('methylationFirebrowse.json',A)



print ('Mutation FIREBROWSE')
Mut2 = firebrowse.Analyses().MutationMAF(format="csv",cohort="COADREAD",tool="MutSig2CV",column="all")
createFile('mutationFirebrowse.csv',Mut2)

print ('Type of Mutations ')
df = pd.read_csv('mutationFirebrowse.csv')
types = df['Variant_Type']

arr = pd.Series.as_matrix(types)
un = np.unique(arr)
print(un)
print(' TOOOOOOOOT ' + str(len(un)))


print ('---------------- miRNA FIREBROWSE ------------------------------------')

miRNA = firebrowse.Samples().miRSeq(cohort='COADREAD',format='csv',mir='tutu')
createFile('miRNAFirebrowse.csv',miRNA)

print (' ----------------- CLINICAL FIREBROWSE ------------------------')
S = firebrowse.Samples().Clinical_FH(cohort="COADREAD",format="csv")
createFile('clinicalFirebrowse.csv',S)

'''
#Demain Question !!!!!
print ('----------------- mRNA FIREBROWSE ---------------------------------------')

mrna = firebrowse.Samples().mRNASeq(cohort='COADREAD',gene='TP53',format='csv',protocol='RPKM')
createFile('mRNAFirebrowseRPKM.csv',mrna)
'''
conn = requests.get(BIGPRX + "getTypesOfCancer")

lines = filterAnswer('col',conn.text)
createFileFromArray('typesOfCancer.csv',lines)

conn = requests.get(BIGPRX + 'getGeneticProfiles&cancer_study_id=coadread_tcga_pub')
createFile('geneticProfiles.csv',conn.text)

print ('--------------- CANCER STUDIES ------------------' )       
conn = requests.get(BIGPRX + "getCancerStudies")
lines = filterAnswer('colorec',conn.text)
createFileFromArray('cancerStudies.csv',lines)

print ('-------------- METHYLATION BIOPORTAL ----------------')
conn = requests.get(BIGPRX + "    "+
"&case_set_id=coadread_tcga_pub_methylation_hm27&"+
"genetic_profile_id=coadread_tcga_pub_methylation_hm27&gene_list=TP53")
createFile('methylationBioPortal.csv',conn.text)

print ('--------- CASE LIST -----------------------')
conn = requests.get(BIGPRX + "getCaseLists&cancer_study_id=coadread_tcga_pub")
createFile('cases.csv',conn.text)
#filterAnswer('hyperme',conn.text)





print ('----------------- mRNA BIOPORTAL --------------------------------------- ')

conn = requests.get(BIGPRX + "getProfileData"+
"&case_set_id=coadread_tcga_pub_rna_seq_mrna&"+
"genetic_profile_id=coadread_tcga_pub_rna_seq_mrna&gene_list=TP53")


createFile('mRNA-BioPortal-KPM.csv',conn.text)

print (' ---------------- miRNA BIOPORTAL -------------------------------------')

microRNA = requests.get(BIGPRX + 'getProfileData&case_set_id=coadread_tcga_pub_microrna&' +
"genetic_profile_id=coadread_tcga_pub_mirna_median_Zscores&gene_list=ARGLU1")

createFile ('miRNAbioportal.csv', microRNA.text)

print ('---------------- CLINICAL BIOPORTAL ---------------------------')

conn = requests.get(BIGPRX + 'getClinicalData&case_set_id=coadread_tcga_pub_all')
createFile('clinicalDataBioPortal.csv',conn.text)
'''
print ('----------------- MUTATION BIOPORTAL --------------------------')

conn = requests.get(BIGPRX + 'getMutationData&case_set_id=coadread_tcga_pub_3way_complete&'+
                    'genetic_profile_id=coadread_tcga_pub_mutations&gene_list=TP53' )
createFile('mutationsBIOPORTAL-CHANGEACTION.csv',conn.text)

GDC_PRX = 'https://gdc-api.nci.nih.gov/'
'''

print ('------------------- GDC RNASEQ colon --------------------------------')
#
# You can use facets to get the possible values that a field can take.
# 'https://gdc-api.nci.nih.gov/cases?facets=files.experimental_strategy&from=1&size=0&pretty=true')
#
cases_endpt = 'https://gdc-api.nci.nih.gov/cases'

filt = {
                    "op":"and",
                    "content":[
                    {
                        "op":"=",
                        "content":{
                                "field": "files.experimental_strategy",
                                "value": "RNA-Seq",
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Colon Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    }
                    ]
            }

params = {'filters':json.dumps(filt),'format':'json','pretty':'true','size':'20'}
# requests URL-encodes automatically
print ('before')
response = requests.get(cases_endpt, params = params)
createFile('FILESMRNAcolon-GDC.json',response.text)


### FILESRMNA is a UUID file list, you have to manually (until now) copy
# paste the file of your choice in :
#
# requests.get('https://gdc-api.nci.nih.gov/data/3a77c136-9ef3-4c0d-baae-6fb98f33d1cd')
#   
#files.experimental_strategy: "miRNA-Seq" , "RNA-Seq" , "Methylation Array"
#project.disease_type="Colon Adenocarcinoma"

conn = requests.get('https://gdc-api.nci.nih.gov/data/3a77c136-9ef3-4c0d-baae-6fb98f33d1cd')

createFile('MRNA-GDCcolon.csv',conn.content)



print ('------------------- GDC RNASEQ rectum --------------------------------')
#
# You can use facets to get the possible values that a field can take.
# 'https://gdc-api.nci.nih.gov/cases?facets=files.experimental_strategy&from=1&size=0&pretty=true')
#
cases_endpt = 'https://gdc-api.nci.nih.gov/cases'

filt = {
                    "op":"and",
                    "content":[
                    {
                        "op":"=",
                        "content":{
                                "field": "files.experimental_strategy",
                                "value": "RNA-Seq",
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Rectum Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    }
                    ]
            }

params = {'filters':json.dumps(filt),'format':'json','pretty':'true','size':'30'}
# requests URL-encodes automatically
print ('before rectum')
response = requests.get(cases_endpt, params = params)
createFile('FILESMRNArectum-GDC.csv',response.text)

### FILESRMNA is a UUID file list, you have to manually (until now) copy
# paste the file of your choice in :
#
# requests.get('https://gdc-api.nci.nih.gov/data/3a77c136-9ef3-4c0d-baae-6fb98f33d1cd')
#   
#files.experimental_strategy: "miRNA-Seq" , "RNA-Seq" , "Methylation Array"
#project.disease_type="Colon Adenocarcinoma"

conn = requests.get('https://gdc-api.nci.nih.gov/data/5807163d-e287-4acf-bb53-35fcc92f5d67')

createFile('MRNA-GDCRectum.csv',conn.content)

print ('-------------------------- GDC miRNA Seq Colon----------------------')

filtMIRNA = {
                    "op":"and",
                    "content":[
                    {
                        "op":"=",
                        "content":{
                                "field": "files.experimental_strategy",
                                "value": "miRNA-Seq",
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Colon Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    }
                    ]
            }

params = {'filters':json.dumps(filtMIRNA),'format':'json','fields': 'file_id','pretty':'true'}
# requests URL-encodes automatically
response = requests.get(cases_endpt, params = params)
createFile('FILESMIRNA-GDC-Colon.csv',response.text)
conn = requests.get('https://gdc-api.nci.nih.gov/data/9593a218-75c6-4495-a48e-ee4b64be4dd3')
createFile('MIRNA-GDC-Colon.csv',conn.content)


print ('-------------------------- GDC miRNA Seq Rectum----------------------')

filtMIRNA = {
                    "op":"and",
                    "content":[
                    {
                        "op":"=",
                        "content":{
                                "field": "files.experimental_strategy",
                                "value": "miRNA-Seq",
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Rectum Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    }
                    ]
            }

params = {'filters':json.dumps(filtMIRNA),'format':'json','fields': 'file_id','pretty':'true','size':'30'}
# requests URL-encodes automatically
response = requests.get(cases_endpt, params = params)
createFile('FILESMIRNA-GDC-Rectum.csv',response.text)
conn = requests.get('https://gdc-api.nci.nih.gov/data/9593a218-75c6-4495-a48e-ee4b64be4dd3')
createFile('MIRNA-GDC-Rectum.csv',conn.content)
 


print ('-------------------------- GDC Colon Methylation Array ----------------------')

filtMeth = {
                    "op":"and",
                    "content":[
                    {
                        "op":"=",
                        "content":{
                                "field": "files.experimental_strategy",
                                "value": "Methylation Array",
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Colon Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    }
                    ]
            }

params = {'filters':json.dumps(filtMeth),'format':'json','fields': 'file_id','pretty':'true'}
# requests URL-encodes automatically
response = requests.get(cases_endpt, params = params)
createFile('FILESMETH-Colon-GDC.csv',response.text)

print ('---- Downloading File colon methylation GDC--------------------')

conn = requests.get('https://gdc-api.nci.nih.gov/data/de7de1aa-d590-4502-b31f-05c0ee6de68a')
createFile('Methylation-GDC-Colon.csv',conn.content)


print ('-------------------------- GDC Rectum Methylation Array ----------------------')

filtMeth = {
                    "op":"and",
                    "content":[
                    {
                        "op":"=",
                        "content":{
                                "field": "files.experimental_strategy",
                                "value": "Methylation Array",
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Rectum Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    }
                    ]
            }

params = {'filters':json.dumps(filtMeth),'format':'json','fields': 'file_id','pretty':'true'}
# requests URL-encodes automatically
response = requests.get(cases_endpt, params = params)
createFile('FILESMETH-Rectum-GDC.csv',response.text)

print ('---- Downloading File colon methylation GDC--------------------')

conn = requests.get('https://gdc-api.nci.nih.gov/data/10092e3d-81f3-44f7-acbe-f05c5c46e817')
createFile('Methylation-GDC-Rectum.csv',conn.content)






#conn = requests.get('https://gdc-api.nci.nih.gov/data/de7de1aa-d590-4502-b31f-05c0ee6de68a')
#createFile('METH.csv',conn.content)

print ('----- clinical Colon GDC ----------')
filtClinical = {
                    "op":"and",
                    "content":[
                    
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Colon Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    },
                        {
                        "op":"=",
                        "content":{
                                "field":"files.data_category",
                                "value":"Clinical"
                        }
                    }
                    ]
    }

params = {'filters':json.dumps(filtClinical),
          'format':'json',
          'pretty':'true'
          }
# requests URL-encodes automatically
response = requests.get('https://gdc-api.nci.nih.gov/cases', params = params)
createFile('CLINICAL-cdgColon.csv',response.text)

print ('-----------------------Clinical Rectum GDC -----------------------')

filtClinical = {
                    "op":"and",
                    "content":[
                    
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Rectum Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    },{
                        "op":"=",
                        "content":{
                                "field":"files.data_category",
                                "value":"Clinical"
                        }
                    }
                    ]
    }

params = {'filters':json.dumps(filtClinical),
          'format':'json',
         
          'pretty':'true'}
# requests URL-encodes automatically
response = requests.get('https://gdc-api.nci.nih.gov/cases', params = params)
createFile('CLINICAL-cdgRectum.csv',response.text)


print ('------------------ Mutations (SNV) Colon GDC-----------------------')


filtClinical = {
                    "op":"and",
                    "content":[
                    
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Colon Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.data_category",
                                "value":"Simple Nucleotide Variation"
                        }
                    }  
                    ]
    }

params = {'filters':json.dumps(filtClinical),
          'format':'json','pretty':'true'}
# requests URL-encodes automatically
response = requests.get('https://gdc-api.nci.nih.gov/cases', params = params)
createFile('mutations-SNV-ColonGDC.csv',response.text)

print ('------------------ Mutations (Copy Number Variation) Colon GDC-----------------------')


filtClinical = {
                    "op":"and",
                    "content":[
                    
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Colon Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.data_category",
                                "value":"Copy Number Variation"
                        }
                    }  
                    ]
    }

params = {'filters':json.dumps(filtClinical),
          'format':'json','pretty':'true'
          }
# requests URL-encodes automatically
response = requests.get('https://gdc-api.nci.nih.gov/cases', params = params)
createFile('mutations-CNV-ColonGDC.csv',response.text)




print ('------------------ Mutations SNV Rectum GDC-----------------------')


filtClinical = {
                    "op":"and",
                    "content":[
                    
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Rectum Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.data_category",
                                "value":"Simple Nucleotide Variation"
                        }
                    }  
                    ]
    }

params = {'filters':json.dumps(filtClinical),
          'format':'json',
          
          'pretty':'true'}
# requests URL-encodes automatically
response = requests.get('https://gdc-api.nci.nih.gov/cases', params = params)
createFile('mutations-SNV-RectumGDC.csv',response.text)

print ('------------------ Mutations CNV Rectum GDC-----------------------')


filtClinical = {
                    "op":"and",
                    "content":[
                    
                    {
                        "op":"=",
                        "content":{
                                "field":"files.cases.project.disease_type",
                                "value":"Rectum Adenocarcinoma"
                        }
                    }, 
                    {
                        "op":"=",
                        "content":{
                                "field":"files.access",
                                "value":"open"
                        }
                    },
                    {
                        "op":"=",
                        "content":{
                                "field":"files.data_category",
                                "value":"Copy Number Variation"
                        }
                    }  
                    ]
    }

params = {'filters':json.dumps(filtClinical),
          'format':'json',
          'pretty':'true'}
# requests URL-encodes automatically
response = requests.get('https://gdc-api.nci.nih.gov/cases', params = params)
createFile('mutations-CNV-Rectum-GDC.csv',response.text)




print ('----------- POSTPROCESS MUTATIONS FIREBROWSE ---------------')

df = pd.read_csv('data/METH27FROZEN/METH27FROZEN.csv')
patientIds = df.iloc[[0]]
print (patientIds)

substr = patientIds.str[8:12]
print(type(substr))
arr = pd.Series.as_matrix(substr)
un = np.unique(arr)
print(len(un))

'''





#print ('***************** FIREBROWSE RNASEQ FROM ARCHIVES *******')
#rnaSEQ = firebrowse.Archives().StandardData(protocol='RSEM_genes_normalized',
#                            cohort='COADREAD',
#                            data_type='mRNASeq',
#                            date='2016_01_28,2015_11_01,2015_08_21,2015_06_01,2015_04_02,2015_02_04,2014_12_06,2014_10_17,2014_09_02,2014_07_15,2014_05_18,2014_04_16,2014_03_16')
#    tool='Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data')
#createFile('RNAAll_____FirebrowseAllDates.json',rnaSEQ)

'''
print('**********************  FIREBROWSE miRNASew from ARCHIVES       ********************************')

rnaSEQ = firebrowse.Archives().StandardData(cohort='COADREAD',data_type='miRSeq',protocol='miR_gene_expression')
#    tool='Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data')
createFile('miRNAAllGeneExpression_____Firebrowse.csv',rnaSEQ)

'''

print ('***************** FIREBROWSE RNASEQ FROM ARCHIVES *******')
rnaSEQ = firebrowse.Archives().StandardData(tool='miRseq_Mature_Preprocess',
                            cohort='COADREAD',
                            data_type='miRSeq',
                            date='2016_01_28,2015_11_01,2015_08_21,2015_06_01,2015_04_02,2015_02_04,2014_12_06,2014_10_17,2014_09_02,2014_07_15,2014_05_18,2014_04_16,2014_03_16')
#    tool='Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data')
createFile('MIRNAAll_____FirebrowseAllDatesMature.json',rnaSEQ)
