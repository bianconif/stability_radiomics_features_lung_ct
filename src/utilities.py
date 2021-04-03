from collections import OrderedDict
from os import listdir
from os.path import isfile, join, splitext

import sqlite3
from dicom_parser import Image

def metadata_from_dicom_folder(dicom_folder):
    """Basic metadata from DICOM folder
    
    Parameters
    ----------
    dicom_folder : str
        The folder where the DICOM files are stored
        
    Returns
    -------
    metadata : OrderedDict
        The metadata (keys describe the meaning).
    """
    
    #Retrieve the list of files in the folder
    dicom_files = [f for f in listdir(dicom_folder) if isfile(join(dicom_folder, f))]
    full_paths = []
    for dicom_file in dicom_files:
        full_paths.append(join(dicom_folder, dicom_file))
    
    #Read the metadata from the first dicom file in the list
    image = None
    for full_path in full_paths:
        _, ext = splitext(full_path)
        if ext in ['.dcm']:
            image = Image(full_path)
            break
    metadata = OrderedDict(
        {'age' : image.header.get('PatientAge'),
         'gender' : image.header.get('PatientSex'),
         'tube_voltage' : image.header.get('KVP')}
    )
    
    return metadata

class DBDriver():
    """Database interface for storing and retrieving the feature values"""
    
    @staticmethod
    def _mangle_feature_name(feature_name):
        return feature_name.replace("/","_")
    
    @staticmethod
    def _experimental_condition(patient_id, nodule_id, annotation_id, 
                                num_levels, noise_scale):
        """Generate string for SQL query that dientifies one combination
        of patient_id, nodule_id, annotation_id, num_levels and noise_scale"""
        condition = f"patient_id='{patient_id}' AND "+\
                    f"nodule_id={nodule_id} AND "+\
                    f"annotation_id={annotation_id} AND "+\
                    f"num_levels={num_levels} AND "+\
                    f"noise_scale={noise_scale}"    
        return condition
    
    def _get_row(self, patient_id, nodule_id, annotation_id, num_levels,
                 noise_scale):
        """Returns the row matching the given patient_id, nodule_id, 
        annotation_id, num_levels and noise_scale"""
        
        condition = self.__class__._experimental_condition(
            patient_id, nodule_id, annotation_id, num_levels, noise_scale)
        command_str = f"SELECT * FROM features WHERE {condition}"
        cur = self._connection.cursor()
        cur.execute(command_str) 
        rows = cur.fetchall()
        
        if len(rows) > 1:
            raise Exception('Internal error: detected multiple rows for the'
                            'same experimental condition') 
        return rows
   
    def _create_new(self):
        """Generates an empty table"""
        
        self._connection = sqlite3.connect(self._db_file)
        
        #Define the table fields
        command_str = "CREATE TABLE features (patient_id text, "+\
                      "nodule_id integer, annotation_id integer, "+\
                      "num_levels integer, noise_scale real"
        feature_cols = ""
        for feature_name in self._feature_names:
            feature_name_modif = self.__class__._mangle_feature_name(
                feature_name)
            feature_cols = feature_cols + f", {feature_name_modif} real"
        command_str = command_str + feature_cols + ')'
            
        #Create the table and commit the changes
        cur = self._connection.cursor()
        cur.execute(command_str)
        self._connection.commit()  
    
    def read_feature_value(self, patient_id, nodule_id, annotation_id, 
                            num_levels, noise_scale, feature_name):
        """Reads one feature value from the database.
        
        Parameters
        ----------
        patient_id : str 
            The patient id.
        nodule_id : int 
            The nodule id
        annotation_id : int
            The annotation id.
        num_levels : int [> 0] 
            The number of quantisation levels
        noise_scale : float 
            The noise scale.
        feature_name : str
            The name of the feature to retrieve.
        
        Returns
        -------
        feature_value : float
            The feature value. None is returned if the feature is not in the
            database.
        """ 
        
        feature_value = None
        feature_name = self.__class__._mangle_feature_name(feature_name)
        command_str = f"SELECT {feature_name} FROM features WHERE "+\
            self.__class__._experimental_condition(
                patient_id, nodule_id, annotation_id, num_levels, noise_scale)
        cur = self._connection.cursor()
        cur.execute(command_str)        
        rows = cur.fetchall()
        
        #Make sure that only one value is returned and raise an exception
        #otherwise
        if len(rows) > 0:
            length_correct = True
            if len(rows) != 1:
                length_correct = False
            else:
                if len(rows[0]) != 1:
                    length_correct = False
            if not length_correct:
                raise Exception('Internal database error found more than one'
                                ' entry for this feature')
            feature_value = rows[0][0]
        
        return feature_value
        
    def write_feature_value(self, patient_id, nodule_id, annotation_id, 
                            num_levels, noise_scale, feature_name, 
                            feature_value):
        """Writes one feature value into the database
        
        Parameters
        ----------
        patient_id : str 
            The patient id.
        nodule_id : int 
            The nodule id
        annotation_id : int
            The annotation id.
        num_levels : int [> 0] 
            The number of quantisation levels
        noise_scale : float 
            The noise scale.
        feature_name : str
            The name of the feature to retrieve.
        """
        
        #Add a new row if necessary or update a field value if the row already 
        #exists
        row = self._get_row(patient_id, nodule_id, annotation_id, num_levels, 
                            noise_scale)
        feature_name = self.__class__._mangle_feature_name(feature_name)
        command_str = None        
        
        if row:
            #Update field in the corresponding row
            condition = self.__class__._experimental_condition(
                patient_id, nodule_id, annotation_id, num_levels, noise_scale)
            command_str = f"UPDATE features SET {feature_name}={feature_value} "+\
                          f"WHERE {condition}"
        else:    
            #Create a new row
            fields = "INSERT INTO features (patient_id, nodule_id, "+\
                     "annotation_id, num_levels, noise_scale, "+\
                    f"{feature_name}) "   
            values = f"VALUES ('{patient_id}', {nodule_id}, "+\
                     f"{annotation_id}, {num_levels}, {noise_scale}, "+\
                     f"{feature_value})"
            command_str = fields + values
        
        cur = self._connection.cursor()
        cur.execute(command_str)
        self._connection.commit() 
           
    def __init__(self, feature_names, db_file):
        """Opens a connection to the db_file if this exists, otherwise creates
        a new file.
        
        Parameters
        ----------
        feature_names : list of str
            The names of the features to compute. 
        num_levels : list of int
            The number of levels for signal quantisation.
        noise_scales : list of float
            The noise intensity.
        db_file : str
            Path to the database file (.db).
        """
        
        self._feature_names = feature_names
        self._db_file = db_file
        self._connection = None
        
        if isfile(db_file):
            self._connection = sqlite3.connect(self._db_file)
        else:
            self._create_new()
            
        
        pass
    