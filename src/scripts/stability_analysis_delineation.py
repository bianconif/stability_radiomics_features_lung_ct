"""Stability of texture features against lesion delineation"""
import pandas as pd
import pingouin as pg
import pylidc as pl

from functions import grade_stability

#Number of requested observers (different lesion delineations) for each nodule
num_requested_annotations = 4

#Noise scale at which the analysis is performed
noise_scale = 0.0

#Number of discretization levels at which the analysis is performed
num_levels = 256

#Get the feature database
df_features = pd.read_csv('cache/feature_db.csv')

#Store the results of stability analysis here
out_file = 'cache/stability_against_delineation.csv'

#Get the list of patients IDs
patient_ids = set(df_features['patient_id'].to_list())

#Select the nodules with the requested number of annotations
selected_nodules = list()
for patient_id in patient_ids:
    nodule_ids = set(
        df_features[df_features['patient_id'] == patient_id]['nodule_id'].to_list()
    )
    for nodule_id in nodule_ids:
        annotation_query = (df_features['patient_id'] == patient_id) &\
                           (df_features['nodule_id'] == nodule_id) &\
                           (df_features['annotation_id'] != -1)
        annotations = df_features[annotation_query]['annotation_id'].to_list()
        
        if len(set(annotations)) == num_requested_annotations:
            selected_nodules.append(
                {'patient_id' : patient_id, 'nodule_id' : nodule_id}
            )

#Get the list of the features available
available_features = set(df_features['feature_name'].to_list())

#Dataframe to store the ICC
df_delineation_icc = pd.DataFrame(columns = ['feature_name', 
                                             'feature_class',
                                             'ICC1k',
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
        annotation_query =\
            (df_features['feature_name'] == feature_name) &\
            (df_features['patient_id'] == selected_nodule['patient_id']) &\
            (df_features['nodule_id'] == selected_nodule['nodule_id']) &\
            (df_features['annotation_id'] != -1)
        annotations = set(df_features[annotation_query]['annotation_id'].to_list()) 
        
        #Fill the data matrix
        data_row = {'feature_name' : feature_name,
                    'patient_id' : selected_nodule['patient_id'],
                    'nodule_id' : selected_nodule['nodule_id'],
                    'patient_and_nodule_id' : selected_nodule_id,
                    }
        for annotation_id, annotation in enumerate(annotations):
            feature_query =\
                (df_features['feature_name'] == feature_name) &\
                (df_features['patient_id'] == selected_nodule['patient_id']) &\
                (df_features['nodule_id'] == selected_nodule['nodule_id']) &\
                (df_features['num_levels'] == num_levels) &\
                (df_features['noise_scale'] == noise_scale) &\
                (df_features['annotation_id'] == annotation)  
            feature_value = df_features[feature_query]['value'].to_list()[0]
            data_row.update({'feature_value' : feature_value,
                             'annotation_id' : annotation})
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
                   'ICC1k' : icc['ICC'][3]}
    df_delineation_icc = df_delineation_icc.append(results_row, 
                                                   ignore_index = True)
df_delineation_icc.to_csv(out_file, index = False)
            
    


            
            