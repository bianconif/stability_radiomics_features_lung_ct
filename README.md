## This is an

## Disclaimer
The information and content available on this repository are provided with no warranty whatsoever. Any use for scientific or any other purpose is conducted at your own risk and under your own responsibility. The authors are not liable for any damages - including any consequential damages - of any kind that may result from the use of the materials or information available on this repository or of any of the products or services hereon described.

## Description and usage
* `src/scripts/compute_features.py` Use this script to compute the features. At the end of the process the features will be stored stored in the `cache/feature.db` file. The calculation may require from a few minutes to several hours depending on the number of features and combinations of parameters requested. The main paramers are:
  - `features_to_compute` the radiomics features to compute (see `feature_lut` in `src/functions` for accepted values);
  - `CT_window` a tuple of two elements (CTmin, CTmax) used to window the input CT signal. Values less or equal than CTmin are set to CTmin; those greater or equal than CTmax to CTmax;
  - `number_of_levelss` a list of positive integers each representing the number of levels used to resample the CT signal between CTmin and CTmax;
  - `noise_scale` a list of non-negative floats each representing the amplitude of the Gaussian noise as a percentage of the spread (standard deviation) of the input signal.
* `src/scripts/patient_population.py` Use this script to retrieve data about the study population at the scan level (_patient id_, _age_, _gender_, _in-plane pixel spacing_, _slice thickness_ and _slice spacing_) and at the nodule level (_patient id_, _nodule id_, _number of annotations_) and at the annotation level (_annotation id_, _subtetly_, _internal structure_, _calcification_, _sphericity_, _margin_, _lobulation_, _spiculation_, _texture_ and _malignancy_ - see [pylidc](https://pylidc.github.io/) documentation for details about the maing of each parameter). The results will be stored in the `cache/scans_metadata.csv` and `cache/noduless_metadata.csv`, respectively.  

## Dependencies
* [NumPy 1.18.5](https://numpy.org/)
* [Pandas 1.1.3](https://pandas.pydata.org/)
* [pylidc 0.2.2](https://pylidc.github.io/)
* [pynrrd 0.4.2](https://pypi.org/project/pynrrd/)
* [pyradiomics 3.0.1](https://pyradiomics.readthedocs.io/en/latest/)
* [SQLite](https://www.sqlite.org/)

## How to cite this work
* Bianconi, F., Palumbo I., Fravolini M.L. _et al_. __Radiomics analysis of lung lesions on CT: experimental evaluation of the stability of texture features against delineation, intensity discretisation and noise__ (to appear)

## References
1.  Armato III, S.G.; McLennan, G.; Bidaut, L. _et al_. __Data From LIDC-IDRI. The Cancer Imaging Archive__. (2015).  http://doi.org/10.7937/K9/TCIA.2015.LO9QL9SX
1. Armato III, S.G., McLennan, G., Bidaut, L. _et al_.; __The Lung Image Database Consortium (LIDC) and Image Database Resource Initiative (IDRI): A completed reference database of lung nodules on CT scans__ (2011) Medical Physics, 38 (2), pp. 915-931. DOI: https://doi.org/10.1118/1.3528204
1. Clark, K., Vendt, B., Smith, K., Freymann, J., Kirby, J., Koppel, P., Moore, S., Phillips, S., Maffitt, D., Pringle, M., Tarbox, L., Prior, F.; __The cancer imaging archive (TCIA): Maintaining and operating a public information repository__
(2013) Journal of Digital Imaging, 26 (6), pp. 1045-1057. DOI: https://doi.org/10.1007/s10278-013-9622-7
1. Van Griethuysen, J.J.M., Fedorov, A., Parmar, C.,  _et al_.; __Computational radiomics system to decode the radiographic phenotype__ (2017) Cancer Research, 77 (21), pp. e104-e107. DOI: https://doi.org/10.1158/0008-5472.CAN-17-0339
