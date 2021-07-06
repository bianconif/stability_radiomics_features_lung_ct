# Impact of lesion delineation and intensity quantisation on the stability of texture features from lung nodules on CT: a reproducible study

## Description and usage

This repository contains the code required for reproducing the results published in the following paper (in the remainder 'the paper'):
* Bianconi, F., Palumbo I., Fravolini M.L., Palumbo, I., Pascoletti, G., Nuvoli, S., Rondini, M., Spanu, A. and Palumbo, B. [Impact of lesion delineation and intensity quantization on the stability of texture features from lung nodules on CT: a reproducible study](https://www.mdpi.com/2075-4418/11/7/1224), Diagnostics, 11(7). art. no. 1224 (2021)

### Computing the features
* Run the `src/scripts/compute_features.py` to compute the texture features. The results will be stored into the `features` table within the `cache/feature.db` file (use [SQLite](https://www.sqlite.org/index.html) to inspect the content). The calculation may require from a few minutes to several hours depending on the number of features and the combinations of parameters requested.

* The first five columns of the `features` table are organised as follows:
  - `patient_id` (id of the scan/patient);
  - `nodule_id` (id of the lung nodule within the scan);
  - `annotation_id` (id of the manual annotation for the nodule);
  - `num_levelss` (number of quantisation levels used for computing the features - parameter N<sub>g</sub>, see Sec. 2.2 of the paper);
  - `noise_scale` (scale of the Gaussian noise; default is 0.0, i.e., no noise, which is the setting used in the paper)

* Each of the remaining columns is labelled as follows:
  - `[feature_class]_[feature_name]`, where `[feature_class]` indicates the feature class (for instance `firstorder`, `glcm`, etc.) and `[feature_name]` the feature name (for instance `entropy`, `Max`, `Mean`, etc.) After the first five columns there will be as many additional columns as the number of features we request to compute.

* Main parameters of the `src/scripts/compute_features.py` script:
  - `features_to_compute` a list containing the names of the radiomics features to compute (see `feature_lut` in `src/functions` for the list of accepted values; please also refer to [pyradiomics](https://pyradiomics.readthedocs.io/en/latest/) documentation for the corresponding definitions and mathematical formulae);
  - `CT_window` a tuple of two elements (CT<sub>min</sub>, CT<sub>max</sub>) representing the clipping bounds for the CT signal (see Sec. 2.2 of the paper);
  - `number_of_levelss` a list of positive integers each representing the number of levels used for signal quantisation (parameter N<sub>g</sub>; see Sec. 2.2 of the paper).
  - `noise_scale` the scale (standard deviation) of the Gaussian noise (not used in the paper; default is 0.0 - no noise)

### Assessing stability against lesion delineation

* Run the `src/scripts/stability_analysis_delineation.py` to assess the stability of the features against lesion delineation. The results will be stored in the `cache/stability_against_delineation.csv` file. The main parameters of the script are:
  - `num_requested_annotations` limits the analysis to those nodules that have exactly the requested number of annotations (default = 4);
  - `num_levelss` (see above);
  - `noise_scale` (see above).

### Assessing stability against intensity resampling

* Run the `src/scripts/stability_analysis_resampling.py` to assess the stability if the features to intensity resampling. The results will be stored in the `cache/stability_against_resampling.csv` file. The main parameters of the script are:
  - `annotation_id` the id of the annotation selected for the analysis (use -1 for 50% consensus - see Sec. 2.4 of the paper);
  - `num_levelss` (see above);
  - `noise_scale` (see above).

### Retrieving the population metadata

* Run the `src/scripts/patient_population.py` to retrieve the data about the study population at the following levels:
  - scan (_patient id_, _age_, _gender_, _in-plane pixel spacing_, _slice thickness_ and _slice spacing_);
  - nodule (_patient id_, _nodule id_, _number of annotations_) and annotation level (_annotation id_, _subtetly_, _internal structure_, _calcification_, _sphericity_, _margin_, _lobulation_, _spiculation_, _texture_ and _malignancy)

The results will be stored in the `cache/scans_metadata.csv` and `cache/nodules_metadata.csv` files, respectively. Please refer to [pylidc](https://pylidc.github.io/) documentation for details about the meaning of each parameter.  

## Python version and dependencies
Tested on Python 3.8.6. Dependencies:
* [NumPy 1.18.5](https://numpy.org/)
* [Pandas 1.1.3](https://pandas.pydata.org/)
* [pylidc 0.2.2](https://pylidc.github.io/)
* [pynrrd 0.4.2](https://pypi.org/project/pynrrd/)
* [pyradiomics 3.0.1](https://pyradiomics.readthedocs.io/en/latest/)
* [SQLite](https://www.sqlite.org/)



## References
1.  Armato III, S.G.; McLennan, G.; Bidaut, L. _et al_. __Data From LIDC-IDRI. The Cancer Imaging Archive__. (2015).  http://doi.org/10.7937/K9/TCIA.2015.LO9QL9SX
1. Armato III, S.G., McLennan, G., Bidaut, L. _et al_.; __The Lung Image Database Consortium (LIDC) and Image Database Resource Initiative (IDRI): A completed reference database of lung nodules on CT scans__ (2011) Medical Physics, 38 (2), pp. 915-931. DOI: https://doi.org/10.1118/1.3528204
1. Clark, K., Vendt, B., Smith, K., Freymann, J., Kirby, J., Koppel, P., Moore, S., Phillips, S., Maffitt, D., Pringle, M., Tarbox, L., Prior, F.; __The cancer imaging archive (TCIA): Maintaining and operating a public information repository__
(2013) Journal of Digital Imaging, 26 (6), pp. 1045-1057. DOI: https://doi.org/10.1007/s10278-013-9622-7
1. Van Griethuysen, J.J.M., Fedorov, A., Parmar, C.,  _et al_.; __Computational radiomics system to decode the radiographic phenotype__ (2017) Cancer Research, 77 (21), pp. e104-e107. DOI: https://doi.org/10.1158/0008-5472.CAN-17-0339

## Acknowledgements
The authors wish to acknowledge the National Cancer Institute and the Foundation for the National Institutes of Health, and their critical role in the creation of the free publicly available [LIDC/IDRI](http://doi.org/10.7937/K9/TCIA.2015.LO9QL9SX) Database used in this study.

## License and data usage
Released under the [Creative Commons Attribution 3.0 Unported License](https://creativecommons.org/licenses/by/3.0/). Users should also abide the [TCIA Data Usage Policy](https://wiki.cancerimagingarchive.net/display/Public/Data+Usage+Policies+and+Restrictions).  

## Disclaimer
The information and content available on this repository are provided with no warranty whatsoever. Any use for scientific or any other purpose is conducted at your own risk and under your own responsibility. The authors are not liable for any damages - including any consequential damages - of any kind that may result from the use of the materials or information available on this repository or of any of the products or services hereon described.
