"""Characteristics of the patient series"""
from collections import OrderedDict

import pandas as pd
import pylidc as pl

from utilities import metadata_from_dicom_folder

#Store the scans metadata here
cache_folder = 'cache'
scans_metadata_csv = cache_folder + '/scans_metadata.csv'

#Store the nodules metadata here
nodules_metadata_csv = cache_folder + '/nodules_metadata.csv'

#Retrieve all the CT scans
scans = pl.query(pl.Scan)
num_scans = scans.count()

df_scans = pd.DataFrame(columns = ['patient_id', 
                                   'age',
                                   'gender',
                                   'pixel_spacing', 
                                   'slice_thickness',
                                   'slice_spacing'])

df_nodules = None                                   

#Keep track of elements to remove for incomplete metadata (lack of gender and/or age)
to_remove = []

#*****************************************************************
#************** Get the metadata at the scan level ***************
#*****************************************************************
scans_indices = range(num_scans)
for s in scans_indices:
           
    #Get the metadata accessible through the pylidc interface
    scan_metadata = OrderedDict(
        {'patient_id' : scans[s].patient_id,
         'pixel_spacing' : scans[s].pixel_spacing,
         'slice_thickness' : scans[s].slice_thickness,
         'slice_spacing' : scans[s].slice_spacing}
    )   
    print(f'Reading metadata of scan {scans[s].patient_id}')
    
    #Get additional metadata directly from the DICOM files
    dicom_folder = scans[s].get_path_to_dicom_files()
    additional_scan_metadata = metadata_from_dicom_folder(dicom_folder)
    
    #Skip scans with incomplete metadata about gender/age
    if None not in additional_scan_metadata.values():
        scan_metadata.update(additional_scan_metadata)
        if df_scans is None:
            df_scans = pd.DataFrame(columns = list(scan_metadata.keys()))
        df_scans = df_scans.append(scan_metadata, ignore_index = True)
    else:
        to_remove.append(s)
        
df_scans.to_csv(scans_metadata_csv, index = False)    
valid_indices = [x for x in scans_indices if x not in to_remove]
#*****************************************************************
#*****************************************************************
#*****************************************************************

#*****************************************************************
#************ Get the metadata at the nodule level ***************
#*****************************************************************  
#Get the list of the selected CT scans
patient_population = pd.read_csv('cache/scans_metadata.csv')
selected_scans = patient_population['patient_id'].tolist()

#Iterate through the scans
for patient_id in selected_scans:
    
    scan = pl.query(pl.Scan).filter(pl.Scan.patient_id == patient_id).first()

    #Get all the nodules within this scan
    nodules = scan.cluster_annotations(verbose = False)
    
    #Iterate through the nodules
    num_nodules = len(nodules)
    for n in range(num_nodules):
               
        #Current nodule
        nodule = nodules[n]
        
        #Iterate through the annotations for the current nodule and get the 
        #corresponding signals and masks
        annotations = list(range(len(nodule)))
        annotation_record = OrderedDict({'patient_id' : patient_id,
                                         'nodule_id' : n,
                                         'num_annotations' : len(annotations)})
        for idx_ann, ann in enumerate(annotations): 
            
            print(f'Reading metadata of scan {scan.patient_id}, '
                  f'nodule {n}, annotation_id {idx_ann}')
            
            annotation_record.update(OrderedDict(
                {'annotation_id' : idx_ann,
                 'subtlety' : nodule[ann].subtlety,
                 'internal_structure' : nodule[ann].internalStructure,
                 'calcification' : nodule[ann].calcification,
                 'sphericity' : nodule[ann].sphericity,
                 'margin' : nodule[ann].margin,
                 'lobulation' : nodule[ann].lobulation,
                 'spiculation' : nodule[ann].spiculation,
                 'texture' : nodule[ann].texture,
                 'malignancy' : nodule[ann].malignancy}
            ))
            if df_nodules is None:
                column_names = list(annotation_record.keys())
                df_nodules = pd.DataFrame(columns = column_names)
            else:
                df_nodules = df_nodules.append(annotation_record, 
                                               ignore_index = True)

df_nodules.to_csv(nodules_metadata_csv, index = False) 