from utilities import DBDriver

feature_names = ['firstorder/Entropy', 'firstorder/IQR', 'firstorder/Kurtosis']
db_driver = DBDriver(feature_names = feature_names, 
                     db_file = 'cache/features.db')
db_driver.write_feature_value(patient_id = 'AA-00', 
                              nodule_id = 1, 
                              annotation_id = 0, num_levels = 128, 
                              noise_scale = 0.05, 
                              feature_name = 'firstorder/Entropy', 
                              feature_value = 0.01)
db_driver.write_feature_value(patient_id = 'AA-00', 
                              nodule_id = 1, 
                              annotation_id = 0, num_levels = 128, 
                              noise_scale = 0.05, 
                              feature_name = 'firstorder/Kurtosis', 
                              feature_value = 0.02)
db_driver.write_feature_value(patient_id = 'AA-00', 
                              nodule_id = 2, 
                              annotation_id = 0, num_levels = 128, 
                              noise_scale = 0.05, 
                              feature_name = 'firstorder/Entropy', 
                              feature_value = 0.03)
db_driver.write_feature_value(patient_id = 'AB-00', 
                              nodule_id = 2, 
                              annotation_id = 0, num_levels = 64, 
                              noise_scale = 0.05, 
                              feature_name = 'firstorder/IQR', 
                              feature_value = 0.03)
row_1 = db_driver._get_row(patient_id = 'AA-00', nodule_id = 1, 
                           annotation_id = 0, num_levels = 128, 
                           noise_scale = 0.05)
row_2 = db_driver._get_row(patient_id = 'AA-01', nodule_id = 1, 
                           annotation_id = 0, num_levels = 128, 
                           noise_scale = 0.05)
f1 = db_driver.read_feature_value(
    patient_id = 'AA-00', nodule_id = 2, annotation_id = 0, num_levels = 128, 
    noise_scale = 0.05, feature_name = 'firstorder/Entropy')
f2 = db_driver.read_feature_value(
    patient_id = 'AA-00', nodule_id = 2, annotation_id = 0, num_levels = 128, 
    noise_scale = 0.05, feature_name = 'firstorder/Kurtosis')
f3 = db_driver.read_feature_value(
    patient_id = 'AA-00', nodule_id = 2, annotation_id = 0, num_levels = 100, 
    noise_scale = 0.05, feature_name = 'firstorder/Entropy')
a = 0