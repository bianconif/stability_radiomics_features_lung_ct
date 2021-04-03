"""Extract images tiles and ground truth from the LIDC dataset"""
import matplotlib.pyplot as plt
import matplotlib.animation as manim

import numpy as np
import pandas as pd

import pylidc as pl
from pylidc.utils import consensus
from radiomics import featureextractor
import six

import SimpleITK as sitk
import nrrd

from functions import compute_feature, preprocess_signal

#Cache folder where to store the nodule signal (subset of the whole scan volume
#enclosed by the bounding box of the nodule) and the corresponding mask
cache_folder = 'cache'
signal_cache = cache_folder + '/signal.nrrd'
mask_cache = cache_folder + '/mask.nrrd'

#List of the features to compute
features_to_compute = ['firstorder/Mean']

#Radiographic solidity: internal texture (solid, ground glass, or mixed))
#Only consider nodules above the given threshold
texture_threshold = 4.0

#CT window
ct_window = (-1350, 150)

#Number of levels for signal resampling 
num_levels = 255

#Level of Gaussian noise
noise_level = 0.025

#Only consider CT scan with slice thickness in this range (values in mm)
slice_thickness_bounds = [2.0, 4.0]

#Only consider nodules with axial diameter in this range (values in mm)
nodule_size_bounds = [5, 40]



#Query for all CT scans within the slice thickness bounds
scans = pl.query(pl.Scan).filter(pl.Scan.slice_thickness >= slice_thickness_bounds[0],
                                 pl.Scan.slice_thickness <= slice_thickness_bounds[1])

#Get the number of the returned scans
num_scans = scans.count()

#Instantiate the pyradiomics feature extractor
extractor_settings = {'binWidth' : 25, 'imageType' : 'Original'}
extractor = featureextractor.RadiomicsFeatureExtractor(**extractor_settings)

#Enable first-order features
extractor.disableAllFeatures()
extractor.enableFeatureClassByName('firstorder')
print('Enabled features:\n\t', extractor.enabledFeatures)

#Iterate through the scans
for s in range(num_scans):
    
    #Get the patient record from the metadata file
    patient_id = scans[s].patient_id
    
    #Get the CT scan as a voxel model
    voxel_model = scans[s].to_volume()
    
    #Get all the nodules within this scan
    nodules = scans[s].cluster_annotations()
    
    #Iterate through the nodules
    num_nodules = len(nodules)
    for n in range(num_nodules):
               
        #Current nodule
        nodule = nodules[n]
        a = 0
        
        #Iterate through the nodules and get the corresponding signals and masks
        for a, annotation in enumerate(nodule):
            
            print(f'Patient ID: {patient_id}, nodule: {n}, annotation: {a}')
            
            signal = voxel_model[annotation.bbox()]
            mask = annotation.boolean_mask().astype(np.uint8)
            
            #*****************************************************************
            #******************* Test preprocessing **************************
            #*****************************************************************
            preprocessed_signal_no_noise = preprocess_signal(
                signal_in = signal, 
                window = ct_window, 
                num_levels = num_levels)
            preprocessed_signal_with_noise = preprocess_signal(
                signal_in = signal, 
                window = ct_window, 
                num_levels = num_levels,
                noise_level = noise_level
            )            
            #*****************************************************************
            #*****************************************************************
            #*****************************************************************
                       
            #Save the signal and mask into the cache folder
            nrrd.write(signal_cache, signal)  
            nrrd.write(mask_cache, mask) 
                      
            #Read back signal and mask with sitk and make sure they are the same
            reread_signal = sitk.GetArrayFromImage(sitk.ReadImage(signal_cache))
            reread_mask = sitk.GetArrayFromImage(sitk.ReadImage(mask_cache))
            
            #Roll the axes
            reread_signal = np.moveaxis(a = reread_signal, 
                                        source = [0,1,2], 
                                        destination = [2,1,0])
            reread_mask = np.moveaxis(a = reread_mask, 
                                      source = [0,1,2], 
                                      destination = [2,1,0])            
            
            signal_diff = np.sum(signal - reread_signal)
            mask_diff = np.sum(mask - reread_mask)
            
            #Compute the radiomics features
            for feature in features_to_compute:
                feature_value = compute_feature(
                    feature_name = feature, 
                    path_to_image = signal_cache,
                    path_to_mask = mask_cache
                )
                a = 0
            
            #features = extractor.execute(signal_cache, mask_cache)
            #print('Radiomics features')
            #for key, value in six.iteritems(features):
                #print('\t', key, ':', value)            
            #a = 0

        
        ##Compute the average diameter over the available annotations
        #diameter = np.mean([annotation.diameter for annotation in nodule])
        
        ##Compute the average texture type over the available annotationz
        #texture = np.mean([annotation.texture for annotation in nodule])
        
        ##Filter by diameter and texture type
        #if (diameter >= nodule_size_bounds[0]) &\
           #(diameter <= nodule_size_bounds[1]) &\
           #(texture >= texture_threshold):
            
            #print(f'Texture level:{texture}')
                    
            ## Perform a consensus consolidation and 50% agreement level.
            ## We pad the slices to add context for viewing.
            #cmask, cbbox, masks = consensus(nodules[n], clevel=0.5, 
                                            #pad=[(30,30), (30,30), (0,0)])
            
            ##Get the index of the slice where the area of the lesion is maximal
            #k = np.argmax(np.sum(masks[0], axis = (0,1)))
            
            ###Get the corresponding volume slice
            ##vol_slice = voxel_model[cbbox][:,:,k]
                        
            ###Rescale the slice
            ##vol_slice_rescaled = rescale_bounds[1] * (vol_slice - ct_window[0])\
                ##/ (ct_window[1] - ct_window[0]) 
            ##vol_slice_rescaled[vol_slice_rescaled < rescale_bounds[0]] = rescale_bounds[0]
            ##vol_slice_rescaled[vol_slice_rescaled > rescale_bounds[1]] = rescale_bounds[1]
            ##vol_slice_rescaled = np.round(vol_slice_rescaled).astype(np.uint8)            
            
            ###Rescale the mask
            ##mask_rescaled = rescale_bounds[1] * cmask[:,:,k]
            ##mask_rescaled[mask_rescaled < rescale_bounds[0]] = rescale_bounds[0]
            ##mask_rescaled[mask_rescaled > rescale_bounds[1]] = rescale_bounds[1]
            ##mask_rescaled = mask_rescaled.astype(np.uint8)
                          
            