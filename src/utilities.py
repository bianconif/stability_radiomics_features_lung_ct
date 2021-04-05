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
    
    @classmethod
    def generate_from_file(cls, db_file):
        """Opens a connection to an existing db_file if this exists.
        
        Parameters
        ----------
        db_file : str
            Path to the database file (.db).
        """  
        
        if not isfile(db_file):
            raise Exception('Database file not found')
        
        #Open a connection
        connection = sqlite3.connect(db_file) 
        
        #Get the feature names
        command_str = "SELECT * FROM PRAGMA_TABLE_INFO('features')"
        cur = connection.cursor()
        cur.execute(command_str) 
        rows = cur.fetchall()    
        feature_names = list()
        to_exclude = {'patient_id', 'nodule_id', 'annotation_id',
                      'noise_scale', 'num_levels'}
        for row in rows:
            if row[1] not in to_exclude:
                feature_names.append(DBDriver._unmangle_feature_name(row[1]))
           
        #Close the connection to the database
        connection.close()
        
        #Instantiate and return the DBDriver
        return DBDriver(feature_names, db_file)
    
    
    @staticmethod
    def _mangle_feature_name(feature_name):
        return feature_name.replace("/","_")
    
    @staticmethod
    def _unmangle_feature_name(feature_name):
        return feature_name.replace("_","/")    
    
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
    
    def _execute_query(self, command_str):
        cur = self._connection.cursor()
        cur.execute(command_str)        
        return cur.fetchall()        
    
    def _get_row(self, patient_id, nodule_id, annotation_id, num_levels,
                 noise_scale):
        """Returns the row matching the given patient_id, nodule_id, 
        annotation_id, num_levels and noise_scale"""
        
        condition = self.__class__._experimental_condition(
            patient_id, nodule_id, annotation_id, num_levels, noise_scale)
        command_str = f"SELECT * FROM features WHERE {condition}"
        rows = self._execute_query(command_str)
        
        if len(rows) > 1:
            raise Exception('Internal error: detected multiple rows for the'
                            'same experimental condition') 
        return rows
    
    def get_feature_names(self):
        return self._feature_names
   
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
        rows = self._execute_query(command_str)
        
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
    
    def get_feature_values_by_annotation(self, patient_id, nodule_id,
                                         feature_name, num_levels = 256,
                                         noise_scale = 0.0):
        """For a given patient, nodule and feature name returns the feature
        values for each of the annotations available in the database. The
        50% consenus annotation is excluded.
        
        Parameters
        ----------
        patient_id : str 
            The patient id.
        feature_name : str
            The name of the feature to retrieve.
        nodule_id : int 
            The nodule id.
        num_levels : int [> 0] 
            The number of quantisation levels
        noise_scale : float 
            The noise scale.
        
        Returns
        -------
        feature_values : list of float
            The feature values. These are as many as the number of delineations
            for the given nodule.
        """
        command_str = f"SELECT {self._mangle_feature_name(feature_name)} FROM features "+\
                      f"WHERE patient_id = '{patient_id}' "+\
                      f"AND nodule_id = {nodule_id} "+\
                      f"AND num_levels = {num_levels} "+\
                      f"AND noise_scale = {noise_scale} "+\
                      f"AND annotation_id != -1 "+\
                      f"ORDER BY annotation_id"
        rows = self._execute_query(command_str)
        feature_values = [row[0] for row in rows]
        return feature_values
    
    def get_patients_ids(self):
        """Returns the unique list of patients' ids
        
        Returns
        -------
        patients_ids : list of str
            The unique list of patients' ids
        """
        
        command_str = f"SELECT DISTINCT patient_id FROM features"
        rows = self._execute_query(command_str)    
        patients_ids = [row[0] for row in rows]
        return patients_ids
    
    def get_nodule_ids_by_patient(self, patient_id):
        """Returns the unique list of nodules' ids for the given patient.
        
        Parameters
        ----------
        patient_id : str 
            The patient id.
        
        Returns
        -------
        nodules_ids : list of int
            The unique list of nodule ids for the given patient.
        """
        
        command_str = f"SELECT DISTINCT nodule_id FROM features "+\
                      f"WHERE patient_id='{patient_id}'"
        rows = self._execute_query(command_str)    
        nodules_ids = [row[0] for row in rows]
        return nodules_ids
    
    def get_annotation_ids_by_nodule(self, patient_id, nodule_id):
        """Returns the unique list of annotations' ids for the given patient
        and nodule (the 50% consenus annotation is excluded).
        
        Parameters
        ----------
        patient_id : str 
            The patient id.
        nodule_id : int 
            The nodule id.
        
        Returns
        -------
        nodules_ids : list of int
            The unique list of nodule ids for the given patient.
        """
        
        command_str = f"SELECT DISTINCT annotation_id FROM features "+\
                      f"WHERE patient_id='{patient_id}' "+\
                      f"AND nodule_id='{nodule_id}' "+\
                      f"AND annotation_id != -1 "
        rows = self._execute_query(command_str)    
        annotation_ids = [row[0] for row in rows]
        return annotation_ids    
        
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
    