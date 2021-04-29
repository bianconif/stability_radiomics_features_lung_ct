"""Stability of texture features against resampling"""
import numpy as np
import pandas as pd
import pingouin as pg
import pylidc as pl

from functions import grade_stability
from utilities import DBDriver

#Number of discretization levels at which the analysis is performed
num_levelss = [32, 64, 128, 256]

#Annotation id for which the analysis is performed
annotation_id = -1

#Noise scale at which the analysis is performed
noise_scale = 0.0

#Get the feature database
db_file = 'cache/features.db'
db_driver = DBDriver.generate_from_file(db_file)

#Store the results of the stability analysis here
out_file = 'cache/stability_against_resampling.csv'

#Get the list of patients IDs
patient_ids = db_driver.get_patients_ids()

#Get the list of the features available
available_features = db_driver.get_feature_names()

#Dataframe to store the ICC
df_noise_icc = pd.DataFrame(columns = ['feature_name', 'feature_class',
                                       'ICC', 'stability'])

#Iterate through the available features and compute the icc for each of them
for feature_name in available_features:
    
    #Data matrix to compute the ICC. Rows represent nodules, columns the
    #feature values corresponding to each level of noise
    noise_columns = ['feature_name', 'patient_id', 'nodule_id',
                     'patient_and_nodule_id', 'noise_scale_id']
    df_data_matrix = pd.DataFrame(columns = noise_columns)    
    
    for patient_id in patient_ids:
        nodule_ids = db_driver.get_nodule_ids_by_patient(patient_id)
        for nodule_id in nodule_ids:        
        
        #Fill the data matrix
            data_row = {'feature_name' : feature_name,
                        'patient_id' : patient_id,
                        'nodule_id' : nodule_id,
                        'patient_and_nodule_id' : f'{patient_id}--{nodule_id}'
                        }
            
            for num_levels_id, num_levels in enumerate(num_levelss):
                feature_value = db_driver.\
                    get_feature_value_on_consensus_annotation(
                        patient_id, nodule_id, feature_name, num_levels, 
                        noise_scale)
                data_row.update({'feature_value' : feature_value,
                                 'num_levels_id' : f'{num_levels_id:d}'})
                df_data_matrix = df_data_matrix.append(data_row, 
                                                       ignore_index = True)
    
    #Convert patient_and_nodule_id from unique strings to unique ints
    _, _, patient_and_nodule_num_id = np.unique(
        df_data_matrix['patient_and_nodule_id'].tolist(),
        return_index = True, return_inverse = True) 
    df_data_matrix['patient_and_nodule_num_id'] = patient_and_nodule_num_id
        
    #Compute the intra-class correlation coefficient
    icc = pg.intraclass_corr(data = df_data_matrix, 
                             targets = 'patient_and_nodule_num_id', 
                             raters = 'num_levels_id', 
                             ratings = 'feature_value').round(3)
    results_row = {'feature_class' : feature_name.split('/', 1)[0],
                   'feature_name' : feature_name.split('/', 1)[1],
                   'stability' : grade_stability(icc['ICC'][1]),
                   'ICC' : icc['ICC'][1]}
    print(results_row)
    df_noise_icc = df_noise_icc.append(results_row, ignore_index = True)
df_noise_icc.to_csv(out_file, index = False)
            
    


            
            