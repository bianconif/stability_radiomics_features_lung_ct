"""Compute the features"""
import os

import pandas as pd
import pylidc as pl
import PySimpleGUI as sg

from functions import get_feature_values
from utilities import DBDriver


#*******************************************************************************
#**************************** Parameters ***************************************
#*******************************************************************************
#Cache folder where to store the nodule signal (subset of the whole scan volume
#enclosed by the bounding box of the nodule) and the corresponding mask
cache_folder = 'cache'
signal_cache = cache_folder + '/signal.nrrd'
mask_cache = cache_folder + '/mask.nrrd'
feature_db = cache_folder + '/features.db'

#Create the cache folder if it doesn't exist
if not os.path.isdir(cache_folder):
    os.makedirs(name = cache_folder)

#List of the features to compute
first_order_statistics = ['Entropy', 'IQR', 'Kurtosis', 'MAD', 'Max', 'Mean', 
                          'Median', 'Min', 'Range', 'RMAD', 'Skewness', 'Std',
                          'Uniformity']
#first_order_statistics = ['Entropy', 'IQR']                          
first_order_statistics = ['firstorder/' + x for x in first_order_statistics]

#Gray Level Co-occurrence Matrix (GLCM) Features
glcm = ['Acorr', 'JointAvg', 'ClProm', 'ClShade', 'ClTen', 'Contrast', 
        'Correlation', 'DiffAvg', 'DiffEnt', 'DiffVar', 'JointEnergy', 
        'JointEntropy', 'IMC1', 'IMC2', 'MCC', 'IDMN', 'ID', 'IDN', 'InvVar', 
        'IDM', 'MaxProb', 'SumAvg',  'SumEnt', 'SumSquares']
glcm = ['glcm/' + x for x in glcm]
   
#Gray Level Dependence Matrix (GLDM) Features
gldm = ['SDE', 'LDE', 'GLN', 'DN', 'DNN', 'GLV', 'DV', 'DE', 'LGLE', 'HGLE',
        'SDLGLE', 'SDHGLE', 'LDLGLE', 'LDHGLE']
gldm = ['gldm/' + x for x in gldm]

#Gray Level Run Length Matrix (GLRLM) Features
glrlm = ['SRE', 'LRE', 'GLN', 'GLNN', 'RLN', 'RLNN', 'RP', 'GLV', 'RV', 'RE',
         'LGLRE', 'HGLRE', 'SRLGLE', 'SRHGLE', 'LRLGLE', 'LRHGLE']
glrlm = ['glrlm/' + x for x in glrlm]

#Gray Level Size Zone Matrix (GLSZM) Features
glszm = ['SAE', 'LAE', 'GLN', 'GLNN', 'SZN', 'SZNN', 'ZP', 'GLV', 'ZV', 'ZE',
         'LGLZE', 'HGLZE', 'SALGLE', 'SAHGLE', 'LALGLE', 'LAHGLE']
glszm = ['glszm/' + x for x in glszm]

#Neighbouring Gray Tone Difference Matrix (NGTDM) Features
ngtdm = ['Coarseness', 'Contrast', 'Busyness', 'Complexity', 'Strength']
ngtdm = ['ngtdm/' + x for x in ngtdm]


#Put all the features together
features_to_compute = first_order_statistics + glcm + gldm + glrlm + glszm +\
    ngtdm

#CT window
ct_window = (-583, 137)

#Number of levels for signal resampling 
num_levelss = [32, 64, 128, 256]

#Level of Gaussian noise
noise_scales = [0.0]

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

#Get the list of the selected CT scans
patient_population = pd.read_csv('cache/scans_metadata.csv')
selected_scans = patient_population['patient_id'].tolist()

#*******************************************************************************
#**************************** Progress window **********************************
#*******************************************************************************

#sg.theme('Dark Red')

BAR_MAX = 100

# layout the Window
layout = [[sg.Text('Patient:'), sg.Text(size = (15,1), key='-pid-')],
          [sg.ProgressBar(BAR_MAX, 
                          orientation='h', 
                          size=(20,20), 
                          key='-patient-progress-')],
          [sg.Text('Nodule:'), sg.Text(size = (3,1), key='-nid-')],
          [sg.ProgressBar(BAR_MAX, 
                          orientation='h', 
                          size=(20,20), 
                          key='-nodule-progress-')],  
          [sg.Text('Annotation:'), sg.Text(size = (3,1), key='-aid-')],
          [sg.ProgressBar(BAR_MAX, 
                          orientation='h', 
                          size=(20,20), 
                          key='-annotation-progress-')],           
          [sg.Text('Noise level: '), sg.Text(size = (5,1), key='-noise-')],
          [sg.Text('Resampling levels: '), sg.Text(size = (5,1), key='-numlev-')],
          [sg.Cancel()]]

# create the Window
window = sg.Window('Custom Progress Meter', layout)
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************

#Create the databse driver
db_driver = DBDriver(feature_names = features_to_compute, 
                     db_file = feature_db)

#Iterate through the scans
for num_patient, patient_id in enumerate(selected_scans):
        
    scan = pl.query(pl.Scan).filter(pl.Scan.patient_id == patient_id).first()
        
    #Get the CT scan as a voxel model
    voxel_model = scan.to_volume()
    
    #Get all the nodules within this scan
    nodules = scan.cluster_annotations(verbose = False)
    
    #Iterate through the nodules
    num_nodules = len(nodules)
    for n in range(num_nodules):
        
        print(f'Patient {num_patient + 1} of {len(selected_scans)}; '
              f'nodule {n} of {num_nodules}')
               
        #Current nodule
        nodule = nodules[n]
        
        #Iterate through the annotations for the current nodule and get the 
        #corresponding signals and masks
        annotations = list(range(len(nodule)))
        annotations.append(-1)              #Add the 50% consensus annotation
        for aid, ann in enumerate(annotations):
                        
            #Get the  value for the requested number of levels and
            #noise scale
            for num_levels in num_levelss:
                for noise_scale in noise_scales:
                    _ = get_feature_values(
                        feature_names = features_to_compute,
                        patient_id = patient_id, 
                        nodule_id = n, 
                        annotation_id = ann, 
                        db_driver = db_driver, 
                        window = ct_window, 
                        num_levels = num_levels, 
                        noise_scale = noise_scale,
                        path_to_image = signal_cache,
                        path_to_mask = mask_cache,
                        verbose=True
                    )  
                    
                    #Update the progress bar
                    event, values = window.read(timeout=10)
                    if event == 'Cancel' or event == sg.WIN_CLOSED:
                        break
                        # update bar with loop value +1 so that bar eventually reaches the maximum
                    window['-patient-progress-'].update(100 * (num_patient + 1)/len(selected_scans))
                    window['-pid-'].update(f'{patient_id}')
                    window['-nodule-progress-'].update(100 * (n + 1)/num_nodules) 
                    window['-nid-'].update(f'{n}')
                    window['-annotation-progress-'].update(100 * (aid + 1)/len(annotations)) 
                    window['-aid-'].update(f'{ann}')                    
                    window['-noise-'].update("{:.1f}%".format(noise_scale)) 
                    window['-numlev-'].update(f'{num_levels}')
                    
window.close()