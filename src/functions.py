import os
import warnings

from itertools import permutations

import numpy as np
import nrrd
import pandas as pd
import pylidc as pl
from pylidc.utils import consensus
from radiomics import featureextractor

feature_lut = {'firstorder/Energy' : {'firstorder' : ['Energy']},
               'firstorder/Entropy' : {'firstorder' : ['Entropy']},
               'firstorder/IQR' : {'firstorder' : ['InterquartileRange']},
               'firstorder/Kurtosis' : {'firstorder' : ['Kurtosis']},               
               'firstorder/MAD' : {'firstorder' : ['MeanAbsoluteDeviation']},
               'firstorder/Mean' : {'firstorder' : ['Mean']},
               'firstorder/Median' : {'firstorder' : ['Median']},
               'firstorder/Max' : {'firstorder' : ['Maximum']},
               'firstorder/Min' : {'firstorder' : ['Minimum']},
               'firstorder/Range' : {'firstorder' : ['Range']},
               'firstorder/RMAD' : {'firstorder' : ['RobustMeanAbsoluteDeviation']},
               'firstorder/Std' : {'firstorder' : ['StandardDeviation']},
               'firstorder/Skewness' : {'firstorder' : ['Skewness']},
               'firstorder/Uniformity' : {'firstorder' : ['Uniformity']},
               'glcm/Acorr' : {'glcm' : ['Autocorrelation']},
               'glcm/JointAvg' : {'glcm' : ['JointAverage']},
               'glcm/ClProm' : {'glcm' : ['ClusterProminence']},
               'glcm/ClShade' : {'glcm' : ['ClusterShade']},
               'glcm/ClTen' : {'glcm' : ['ClusterTendency']},
               'glcm/Contrast' : {'glcm' : ['Contrast']},
               'glcm/Correlation' : {'glcm' : ['Correlation']},
               'glcm/DiffAvg' : {'glcm' : ['DifferenceAverage']},
               'glcm/DiffEnt' : {'glcm' : ['DifferenceEntropy']},
               'glcm/DiffVar' : {'glcm' : ['DifferenceVariance']},
               'glcm/JointEnergy' : {'glcm' : ['JointEnergy']},
               'glcm/JointEntropy' : {'glcm' : ['JointEntropy']},
               'glcm/IMC1' : {'glcm' : ['Imc1']},
               'glcm/IMC2' : {'glcm' : ['Imc2']},   
               'glcm/MCC' : {'glcm' : ['MCC']},
               'glcm/IDMN' : {'glcm' : ['Idmn']},
               'glcm/ID' : {'glcm' : ['Id']},
               'glcm/IDN' : {'glcm' : ['Idn']},
               'glcm/InvVar' : {'glcm' : ['InverseVariance']},
               'glcm/IDM' : {'glcm' : ['Idm']},
               'glcm/MaxProb' : {'glcm' : ['MaximumProbability']},
               'glcm/SumAvg' : {'glcm' : ['SumAverage']},
               'glcm/SumEnt' : {'glcm' : ['SumEntropy']},
               'glcm/SumSquares' : {'glcm' : ['SumSquares']},
               'gldm/SDE' : {'gldm' : ['SmallDependenceEmphasis']},
               'gldm/LDE' : {'gldm' : ['LargeDependenceEmphasis']},
               'gldm/GLN' : {'gldm' : ['GrayLevelNonUniformity']},
               'gldm/DN' : {'gldm' : ['DependenceNonUniformity']},
               'gldm/DNN' : {'gldm' : ['DependenceNonUniformityNormalized']},
               'gldm/GLV' : {'gldm' : ['GrayLevelVariance']},
               'gldm/DV' : {'gldm' : ['DependenceVariance']},
               'gldm/DE' : {'gldm' : ['DependenceEntropy']},
               'gldm/LGLE' : {'gldm' : ['LowGrayLevelEmphasis']},
               'gldm/HGLE' : {'gldm' : ['HighGrayLevelEmphasis']},
               'gldm/SDLGLE' : {'gldm' : ['SmallDependenceLowGrayLevelEmphasis']},
               'gldm/SDHGLE' : {'gldm' : ['SmallDependenceHighGrayLevelEmphasis']},
               'gldm/LDLGLE' : {'gldm' : ['LargeDependenceLowGrayLevelEmphasis']},
               'gldm/LDHGLE' : {'gldm' : ['LargeDependenceHighGrayLevelEmphasis']},
               'glrlm/SRE' : {'glrlm' : ['ShortRunEmphasis']},
               'glrlm/LRE' : {'glrlm' : ['LongRunEmphasis']},
               'glrlm/GLN' : {'glrlm' : ['GrayLevelNonUniformity']},
               'glrlm/GLNN' : {'glrlm' : ['GrayLevelNonUniformityNormalized']},
               'glrlm/RLN' : {'glrlm' : ['RunLengthNonUniformity']},
               'glrlm/RLNN' : {'glrlm' : ['RunLengthNonUniformityNormalized']},
               'glrlm/RP' : {'glrlm' : ['RunPercentage']},
               'glrlm/GLV' : {'glrlm' : ['GrayLevelVariance']},
               'glrlm/RV' : {'glrlm' : ['RunVariance']},
               'glrlm/RE' : {'glrlm' : ['RunEntropy']},
               'glrlm/LGLRE' : {'glrlm' : ['LowGrayLevelRunEmphasis']},
               'glrlm/HGLRE' : {'glrlm' : ['HighGrayLevelRunEmphasis']},
               'glrlm/SRLGLE' : {'glrlm' : ['ShortRunLowGrayLevelEmphasis']},
               'glrlm/SRHGLE' : {'glrlm' : ['ShortRunHighGrayLevelEmphasis']},
               'glrlm/LRLGLE' : {'glrlm' : ['LongRunLowGrayLevelEmphasis']},
               'glrlm/LRHGLE' : {'glrlm' : ['LongRunHighGrayLevelEmphasis']},
               'glszm/SAE' : {'glszm' : ['SmallAreaEmphasis']},
               'glszm/LAE' : {'glszm' : ['LargeAreaEmphasis']},
               'glszm/GLN' : {'glszm' : ['GrayLevelNonUniformity']},
               'glszm/GLNN' : {'glszm' : ['GrayLevelNonUniformityNormalized']},
               'glszm/SZN' : {'glszm' : ['SizeZoneNonUniformity']},
               'glszm/SZNN' : {'glszm' : ['SizeZoneNonUniformityNormalized']},
               'glszm/ZP' : {'glszm' : ['ZonePercentage']},
               'glszm/GLV' : {'glszm' : ['GrayLevelVariance']},
               'glszm/ZV' : {'glszm' : ['ZoneVariance']},
               'glszm/ZE' : {'glszm' : ['ZoneEntropy']},
               'glszm/LGLZE' : {'glszm' : ['LowGrayLevelZoneEmphasis']},
               'glszm/HGLZE' : {'glszm' : ['HighGrayLevelZoneEmphasis']},
               'glszm/SALGLE' : {'glszm' : ['SmallAreaLowGrayLevelEmphasis']},
               'glszm/SAHGLE' : {'glszm' : ['SmallAreaHighGrayLevelEmphasis']},
               'glszm/LALGLE' : {'glszm' : ['LargeAreaLowGrayLevelEmphasis']},
               'glszm/LAHGLE' : {'glszm' : ['LargeAreaHighGrayLevelEmphasis']},
               'ngtdm/Coarseness' : {'ngtdm' : ['Coarseness']},
               'ngtdm/Contrast' : {'ngtdm' : ['Contrast']},
               'ngtdm/Busyness' : {'ngtdm' : ['Busyness']},
               'ngtdm/Complexity' : {'ngtdm' : ['Complexity']},
               'ngtdm/Strength' : {'ngtdm' : ['Strength']},
               'shape/MaxAxialDiameter' : {'shape' : ['Maximum2DDiameterSlice']}
               }

def preprocess_signal(signal_in, window = (-1350, 150), num_levels = 256,
                      **kwargs):
    """CT data preprocessing
    
    Parameters
    ----------
    signal_in : a 3D nparray of int or float 
        The input CT data. May represent a whole scan or a part of it.
    window : a list or tuple of float (lower_bound, upper_bound)
        The window bounds in Hounsfield Units.
    num_levels : int (> 1)
        The number of levels used signal image quantisation (resampling).
    noise_scale : float (> 0.0, optional)
        Scale of the Gaussian noise to be added to the original signal. The value
        indicates the spread (standard deviation) of the noise as a percentage of 
        the spread of the input signal. For instance, use noise_scale = 2.5 to
        add Gaussian noise sampled from a normal distribution with spread =
        0.025 that of the original signal.
     
    Returns
    -------
    signal_out : a 3D array of int16 (same size as image_in)
    """
    
    #Convert the input signal to float
    signal_in = signal_in.astype(np.float)
    

    #Add Gaussian noise if required
    if ('noise_scale' in kwargs.keys()) & (kwargs['noise_scale'] > 0.0):
        noise_scale = kwargs['noise_scale']/100
        
        #Normalise the input signal to zero mean and unit variance
        mean = np.mean(signal_in.flatten())
        std = np.std(signal_in.flatten())
        zero_mean_unit_var = (signal_in - mean)/std
        
        #Generate and add the noise
        noise = np.random.normal(loc = 0.0, scale = noise_scale, 
                                 size = signal_in.shape)
        with_noise = zero_mean_unit_var + noise
        
        #Convert back to the original units
        signal_in = with_noise*std + mean

            
    #Normalise the signal to [0,1] according to the given window
    normalised_image = (signal_in - window[0])/(window[1] - window[0])
    normalised_image[normalised_image < 0.0] = 0.0
    normalised_image[normalised_image > 1.0] = 1.0
    
    #Resample the signal to the number of levels required
    signal_out = np.round(normalised_image*(num_levels - 1))/((num_levels - 1))
    
    #Covert back to the original units
    signal_out = (window[1] - window[0])*signal_out + window[0]
    
    return signal_out
    
def get_feature_values(feature_names, patient_id, nodule_id, annotation_id, 
                       db_driver, window, num_levels, noise_scale, path_to_image, 
                       path_to_mask, verbose=False):
    """Value of a set of radiomic features for a given patient and nodule id. 
    The function parses the csv_cache first to check if all the requested 
    feature values are already there; if so reads the values and returns them, 
    otherwise computes the missing values, updates the cache and returns all
    the requested values.
    
    Parameters
    ----------
    feature_names : list of str
        The names of the features to be computed. Possible values are the keys
        in feature_lut dict.
    patient_id : str
        The patient id.
    nodule_id : int
        The nodule id.
    annotation_id : int
        The annotation id for the given nodule. Each annotation corresponds to
        a different observer. Use -1 for 50% consensus annotation
    db_driver : DBDriver
        Instance of a DBDriver which manages feature caching.
    window : a list or tuple of float (lower_bound, upper_bound)
        The window bounds in Hounsfield Units.
    num_levels : int (> 1)
        The number of levels used signal image quantisation (resampling).
    noise_scale : float (> 0.0)
        Scale of the Gaussian noise to be added to the original signal. The value
        indicates the spread (standard deviation) of the noise as a percentage of 
        the spread of the input signal. For instance, use noise_scale = 2.5 to
        add Gaussian noise sampled from a normal distribution with spread =
        0.025 that of the original signal.
    path_to_image : str
        Path to the temporary file where the signal is to be stored. Needs to 
        be an .nrrd file.
    path_to_mask : str
        Path to the temporary file where the mask is to be stored. Needs to 
        be an .nrrd file.
    verbose : bool
        Print details about the features being computed.
    
    Returns
    -------
    feature_values : list of float
        The values of the requested features
    """
    
    feature_names_and_values = dict()
    
    #Read from the database the values of the features that have already been 
    #computed
    for feature_name in feature_names:
        feature_value = db_driver.read_feature_value(
            patient_id, nodule_id, annotation_id, num_levels, noise_scale, 
            feature_name)

        if feature_value:
            feature_names_and_values.update({feature_name : feature_value})
    
    if verbose:
        if len(feature_names_and_values) > 0:
            record_str = f'patient_id : {patient_id}, nodule_id : {nodule_id}, '+\
                f'annotation_id : {annotation_id}, num_levels : {num_levels} '+\
                f'noise_scale : {noise_scale}\n'+\
                f'feature_names : {list(feature_names_and_values.keys())}'
            print(f'Retrived {record_str}')             
            
    
    #Determine which features still needs to be computed
    names_of_features_to_compute = set(feature_names).\
        difference(set(feature_names_and_values.keys()))
    
    #Compute the features that are not in the cache    
    if len(names_of_features_to_compute) > 0:
        if verbose:
            record_str = f'patient_id : {patient_id}, nodule_id : {nodule_id}, '+\
                f'annotation_id : {annotation_id}, num_levels : {num_levels} '+\
                f'noise_scale : {noise_scale}\n'+\
                f'feature_names : {names_of_features_to_compute}'        
            print(f'Computing {record_str}')
    
        #Get the scan corresponding to the given patient_id
        scan = pl.query(pl.Scan).filter(pl.Scan.patient_id == patient_id).first()
    
        #Get the CT scan as a voxel model
        voxel_model = scan.to_volume()    
    
        #Get the requested nodule within this scan
        nodules = scan.cluster_annotations() 
        nodule = nodules[nodule_id]
    
        #Get the requested annotation (nodule mask) and the corresponding signal
        if annotation_id == -1:
            #Get the 50% consensus annotation
            mask, bbox, _ = consensus(nodule, clevel=0.5)        
        else:
            annotation = nodule[annotation_id]
            bbox = annotation.bbox()
            mask = annotation.boolean_mask() 
        mask = mask.astype(np.uint8)
        signal = voxel_model[bbox]
        
        #Preprocess the signal
        signal = preprocess_signal(signal_in = signal, 
                                   window = window, 
                                   num_levels = num_levels,
                                   noise_scale = noise_scale)        
    
        #Store the signal and mask as temporary files
        nrrd.write(path_to_image, signal)  
        nrrd.write(path_to_mask, mask)     
    
        #Compute the feature values and update the database
        values_of_features_to_compute = compute_feature_values(
            names_of_features_to_compute, path_to_image, path_to_mask,
            bin_width = (window[1] - window[0])/num_levels)
        for f, name_of_features_to_compute in enumerate(names_of_features_to_compute):
            feature_names_and_values.update({name_of_features_to_compute :
                                             values_of_features_to_compute[f]})
            db_driver.write_feature_value(patient_id, nodule_id, annotation_id, 
                                          num_levels, noise_scale, 
                                          name_of_features_to_compute, 
                                          values_of_features_to_compute[f])             
    
    #Arrange the requested features in the correct order
    feature_values = list()
    for feature_name in feature_names:
        feature_values.append(feature_names_and_values[feature_name])    
    
    return feature_values
    
    
def compute_feature_values(feature_names, path_to_image, path_to_mask, 
                           bin_width = 1):
    """Compute a set of radiomic features
    
    Parameters
    ----------
    feature_names : list of str
        The names of the feature to compute. Possible values are tke keys of
        feature_lut dict.
    path_to_image : str
        Path to the source image (signal). Needs to be a .nrrd file.
    path_to_mask : str
        Path to the mask image. Needs to be a .nrrd file.
    bin_width : float
        A positive value representing the size of the bins when making a 
        histogram and for discretization of the image gray level. 
        
    Returns
    -------
    feature_values : list of float
        The values of the requested features 
    """
    
    feature_values = list()
    
    #Make sure the requested feature is in the dictionary of the computable ones
    for feature_name in feature_names:
        if feature_name not in feature_lut.keys():
            raise Exception('Feature name not found in the lookup table')
    
    #Instantiate the feature extractor
    settings = {'binWidth' : bin_width}
    extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
    
    #Enable the extraction of the feature requested
    extractor.disableAllFeatures() 
    classes = set([x.split('/', 1)[0] for x in feature_names])
    features_to_enable = {}
    for class_ in classes:
        features_to_enable[class_] = list()
    for feature_name in feature_names:
        class_ = feature_name.split('/', 1)[0]
        features_to_enable[class_].append(feature_lut[feature_name][class_][0])
    extractor.enableFeaturesByName(**features_to_enable)
    
    #Compute the feature requested
    results = extractor.execute(path_to_image, path_to_mask)
    
    #Retrieve the feature value from the results dictionary
    internal_feature_name = str()
    for feature_name in feature_names:
        for key, value in feature_lut[feature_name].items():
            internal_feature_name = key + '_' + value[0]
            for result_name, result_value in results.items():
                tail = result_name.split('_', 1)[1]
                if internal_feature_name == tail:
                    feature_values.append(result_value.tolist())
                    break
    
    #Deinstantiate the extractor explictly 
    del extractor
    
    return feature_values

def grade_stability(avg_smape):
    """Qualitative label for the average symmetric mean absolute percentage 
    error (SMAPE).
    
    Parameters
    ----------
    avg_smape : float (> 0)
        The average absolute relative difference.
    
    Returns
    -------
    qualitative_label : str
        The qualitative label for the given average SMAPE. Can be: 'poor', 
        'moderate', 'good' or 'excellent'
    """
    
    qualitative_label = None
    
    if avg_smape <= 5.0:
        qualitative_label = 'excellent' 
    elif avg_smape <= 10.0:
        qualitative_label = 'good' 
    elif avg_smape <= 20.0:
        qualitative_label = 'moderate'
    else:
        qualitative_label = 'poor'
    
    return qualitative_label

def smape(a, f):
    """Symmetric mean absolute percentage error (SMAPE). Based on the formula
    on page 406 of [1], but with the factor '2' at the denominator removed
    so as to have results bounded in [0,100].
    
    Parameters
    ----------
    a : iterable of numerics
        The actual (target) values.
    f : iterable of numerics
        The forecasted (predicted) values.
    NOTE: a and f must have the same length.
        
    Returns
    -------
    smape_value : float
        The error value.
        
    References
    ----------
    Goodwin, P., Lawton, R. On the asymmetry of the symmetric MAPE (1999) 
    International Journal of Forecasting, 15 (4), pp. 405-408.
    """
    
    smape_value = 0.0
    
    if len(a) != len(f):
        raise Exception('The actual and forecasted values must have the '
                        'same length')    
    
    #Remove zero values to avoid NaN
    nonzero_a_f = np.intersect1d(a.nonzero(), f.nonzero())
    if len(nonzero_a_f) > 0:
        a = a[nonzero_a_f]
        f = f[nonzero_a_f]
        smape_value = 1/len(a) *\
            np.sum(np.abs(f-a) / (np.abs(a) + np.abs(f))*100)
    else:
        warnings.warn(f"SMAPE: no non-zero values in the input arrays, " +\
                      f"returning default value ({smape_value})")
    return smape_value

def avg_smape(values):
    """Average symmetric mean absolute percentage error over an array of
    values. The values ideally represent the results of repeated mesurements
    on the same subject.
    
    Parameters
    ----------
    values : iterable of numerics.
       The input values.
       
    Returns
    -------
    avg_smape : float
       The average SMAPE.
    """
    
    #Get all the pairwise combinations (order matters) and consider the first
    #element as the target and the second as the forecast
    indices = list(range(len(values)))
    combs = permutations(indices, 2)
    targets = list()
    forecasts = list()
    for comb in combs:
        targets.append(values[comb[0]])
        forecasts.append(values[comb[1]])
        
    avg_smape = smape(np.asarray(targets), np.asarray(forecasts))
    return avg_smape
    
 
    
    