"""Stability of texture features against lesion delineation"""
import pandas as pd
import pingouin as pg
import pylidc as pl

from functions import grade_stability
from utilities import DBDriver

#Number of requested observers (different lesion delineations) for each nodule
num_requested_annotations = 4

#Noise scale at which the analysis is performed
noise_scale = 0.0

#Number of discretization levels at which the analysis is performed
num_levels = 256

#Get the feature database
db_file = 'cache/features.db'
db_driver = DBDriver.generate_from_file(db_file)

#Store the results of the stability analysis here
out_file = 'cache/stability_against_delineation.csv'

#Get the list of patients IDs
patient_ids = db_driver.get_patients_ids()

#Select the nodules with the requested number of annotations
selected_nodules = list()
for patient_id in patient_ids:
    nodule_ids = db_driver.get_nodule_ids_by_patient(patient_id)
    for nodule_id in nodule_ids:
        annotation_ids = db_driver.get_annotation_ids_by_nodule(patient_id, 
                                                                nodule_id)        
        if len(annotation_ids) == num_requested_annotations:
            selected_nodules.append(
                {'patient_id' : patient_id, 'nodule_id' : nodule_id}
            )

#Get the list of the features available
available_features = db_driver.get_feature_names()

#Dataframe to store the ICC
df_delineation_icc = pd.DataFrame(columns = ['feature_name', 
                                             'feature_class',
                                             'ICC',
                                             'stability']
                                  )

#Iterate through the available features and compute the icc for each of them
for feature_name in available_features:
    
    #Data matrix to compute the ICC. Rows represent subjects, columns the
    #feature values corresponding to each delineation
    annotation_columns = ['feature_name', 'patient_id', 'nodule_id',
                          'patient_and_nodule_id',
                          'annotation_id']
    df_data_matrix = pd.DataFrame(columns = annotation_columns)    
    
    for selected_nodule_id, selected_nodule in enumerate(selected_nodules):
        patient_id = selected_nodule['patient_id']
        nodule_id = selected_nodule['nodule_id']
        
        #Fill the data matrix
        data_row = {'feature_name' : feature_name,
                    'patient_id' : patient_id,
                    'nodule_id' : nodule_id,
                    'patient_and_nodule_id' : selected_nodule_id,
                    }

        feature_values = db_driver.get_feature_values_by_annotation(
            patient_id, nodule_id, feature_name, num_levels, noise_scale)
        for annotation_id, feature_value in enumerate(feature_values):
            data_row.update({'feature_value' : feature_value,
                             'annotation_id' : annotation_id})
            df_data_matrix = df_data_matrix.append(data_row, 
                                                   ignore_index = True)
    
    #Compute the intra-class correlation coefficient
    icc = pg.intraclass_corr(data = df_data_matrix, 
                             targets = 'patient_and_nodule_id', 
                             raters = 'annotation_id', 
                             ratings = 'feature_value').round(3)
    results_row = {'feature_class' : feature_name.split('/', 1)[0],
                   'feature_name' : feature_name.split('/', 1)[1],
                   'stability' : grade_stability(icc['ICC'][3]),
                   'ICC' : icc['ICC'][3]}
    print(results_row)
    df_delineation_icc = df_delineation_icc.append(results_row, 
                                                   ignore_index = True)
df_delineation_icc.to_csv(out_file, index = False)
            
    


            
            