"""Stability of texture features against resampling"""
import numpy as np
import pandas as pd
import pylidc as pl

from functions import grade_stability, avg_percent_abs_diff
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
#available_features = ['firstorder/Max']

#Dataframe to store the ICC
df_resampling_stability = pd.DataFrame(columns = ['feature_name', 'feature_class',
                                       'ICC', 'stability'])

#Iterate through the available features and compute the icc for each of them
for feature_name in available_features:
    
    abs_relative_variation = list()
      
    for patient_id in patient_ids:
        nodule_ids = db_driver.get_nodule_ids_by_patient(patient_id)
        for nodule_id in nodule_ids:        
        
            #Store the feature values by number of sampling levels
            feature_values_by_sampling_levels = list()
            
            for num_levels_id, num_levels in enumerate(num_levelss):
                feature_value = db_driver.\
                    get_feature_value_on_consensus_annotation(
                        patient_id, nodule_id, feature_name, num_levels, 
                        noise_scale)
                feature_values_by_sampling_levels.append(feature_value)
            
            #Compute the average relative variation by nodule and append it to 
            #the list
            abs_relative_variation_by_nodule = avg_percent_abs_diff(
                feature_values_by_sampling_levels)
            abs_relative_variation.append(abs_relative_variation_by_nodule)
            a = 0
        
    #Compute the average relative variation for the whole population
    avg_abs_relative_variation = np.mean(abs_relative_variation)

    print(f'Feature: {feature_name}; '
          f'Average absolute variation: {avg_abs_relative_variation:.2f}%')
    
    #Generate record for csv output
    results_row = {'feature_class' : feature_name.split('/', 1)[0],
                   'feature_name' : feature_name.split('/', 1)[1],
                   'stability' : grade_stability(avg_abs_relative_variation),
                   'ICC' : avg_abs_relative_variation}    
    
    df_resampling_stability.append(results_row, ignore_index = True)
df_resampling_stability.to_csv(out_file, index = False)
            
    


            
            