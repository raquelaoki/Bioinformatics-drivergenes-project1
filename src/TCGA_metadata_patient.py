# -*- coding: utf-8 -*-
"""
Working with TCGA metadata
Author: Raquel Aoki
Date: 2017/012/27

This code check if a patient is present in 
5 different types of data from TCGA

"""

'''
Each file was download separated
'''
#Load libraries
import json
import pandas as pd
import numpy as np

#files names downloaded
fn = pd.read_csv('TCGA_filesnames.csv')

'''
Vectors to save size and intersection
'''
mutation_size = []
copynumber_size = []
transcriptomic_size = []
methylation_size =[]
clinical_size = []
inter_mct_size = []
inter_mctc_size = []

dem_gender_female = []
dem_gender_male = []
dem_eth_white = []
dem_eth_black = []
dem_etc_not = []
dem_etc_asian = []
dem_etc_native = []
dem_year_min = []
dem_year_mean = []
dem_year_max = []
dem_year_median = []


for r in range(fn.shape[0]):
    # Reading each file
    #mutation
    with open( fn.base1_name[r], 'r') as f:
         data1 = json.load(f)[0]
    
    #copynumber
    with open( fn.base2_name[r], 'r') as f:
         data2 = json.load(f)
    
    #transcriptomic     
    with open( fn.base3_name[r], 'r') as f:
         data3 = json.load(f)
    
    #methylation
    with open(fn.base4_name[r],'r') as f:
        data4 = json.load(f)

    #clinical
    with open( fn.clinical[r], 'r') as f:
         clinical = json.load(f)
    
    '''
    Saving TCGA codes
    '''
    #mutation
    bd1 = data1['analysis']["input_files"]
    mutation_id = []
    for i in range(len(bd1)):
        mutation_id.append(bd1[i]["submitter_id"][:15])
    mutation_id = np.unique(mutation_id)
    
    #copynumber
    copynumber_id = []
    for i in range(len(data2)):
        bd2 = data2[i]["associated_entities"]
        copynumber_id.append(bd2[0]['entity_submitter_id'][:15])
    
    copynumber_id = np.unique(copynumber_id)
        
    #transcroptomic
    transcriptomic_id = []
    for i in range(len(data3)):
        bd3 = data3[i]["associated_entities"]
        transcriptomic_id.append(bd3[0]['entity_submitter_id'][:15])
    
    transcriptomic_id = np.unique(transcriptomic_id)
    
    #methylation
    methylation_id = []
    for i in range(len(data4)):
        bd4 = data4[i]["associated_entities"]
        methylation_id.append(bd4[0]["entity_submitter_id"][:15])
    
    methylation_id = np.unique(methylation_id)
    
    #clinical
    clinical_id = []
    demographic_gender = []
    demographic_year_of_birth = []
    demographic_ethnicity = []
        
    for i in range(len(clinical)):
        bd3 = clinical[i]["associated_entities"]
        clinical_id.append(bd3[0]['entity_submitter_id'][:12])
        bd4 = clinical[i]["cases"]
        bd4 = bd4[0]['demographic']
        demographic_gender.append(bd4['gender'])
        demographic_year_of_birth.append(bd4['year_of_birth'])
        demographic_ethnicity.append(bd4['race'])
        
    clinical_id = np.unique(clinical_id)
    #creating tables and descritive analysis
    demographic_gender = pd.Series(demographic_gender)
    demographic_year_of_birth = pd.Series(demographic_year_of_birth)
    demographic_ethnicity = pd.Series(demographic_ethnicity)
            
    dem_gender_female.append(sum(demographic_gender=='female'))
    dem_gender_male.append(sum(demographic_gender=='male'))
    dem_eth_white.append(sum(demographic_ethnicity=='white'))
    dem_eth_black.append(sum(demographic_ethnicity=='black or african american'))
    dem_etc_not.append(sum(demographic_ethnicity=='not reported'))
    dem_etc_asian.append(sum(demographic_ethnicity=='asian'))
    dem_etc_native.append(sum(demographic_ethnicity=='american indian or alaska native'))
    dem_year_min.append(demographic_year_of_birth.min())
    dem_year_mean.append(demographic_year_of_birth.mean())
    dem_year_max.append(demographic_year_of_birth.max())
    dem_year_median.append(demographic_year_of_birth.median())

    
    '''
    saving vector sizes
    '''
    print(len(mutation_id),len(copynumber_id),len(transcriptomic_id),len(clinical_id))
    mutation_size.append(len(mutation_id))
    copynumber_size.append(len(copynumber_id))
    transcriptomic_size.append(len(transcriptomic_id))
    clinical_size.append(len(clinical_id))
    methylation_size.append(len(methylation_id))
    '''
    intersection
    '''
    aux1 = set(mutation_id).intersection(copynumber_id)
    aux1 = set(aux1).intersection(transcriptomic_id)
    aux1 = set(aux1).intersection(methylation_id)
    inter_mct_size.append(len(aux1))
    
    '''
    edit code to match with clinical info
    '''
    aux1 = list(aux1)
    aux2 = []
    for i in aux1:
        aux2.append(i[:12])
    
    aux2 = set(aux2).intersection(clinical_id)
    inter_mctc_size.append(len(aux2))


fn['mutation_size'] = mutation_size
fn['copy_size'] = copynumber_size
fn['transc_size'] = transcriptomic_size
fn['inter3'] = inter_mct_size
fn['clinical_size2'] = clinical_size
fn['methylation'] = methylation_size
fn['inter4'] = inter_mctc_size

fn['dem_gender_female'] = dem_gender_female
fn['dem_gender_male'] = dem_gender_male
fn["dem_eth_white"] = dem_eth_white
fn['dem_eth_black'] = dem_eth_black
fn['dem_etc_not'] = dem_etc_not
fn['dem_etc_asian'] = dem_etc_asian
fn['dem_etc_native'] = dem_etc_native
fn['dem_year_min'] = dem_year_min 
fn['dem_year_mean'] = dem_year_mean
fn['dem_year_max'] = dem_year_max
fn['dem_year_median'] = dem_year_median

fn['inter5'] = inter_mctc_size

fn.to_csv('TCGA_datasets_sizes2.csv',sep=',')

