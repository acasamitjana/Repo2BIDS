import pdb
from os.path import join, dirname, exists
from os import makedirs
import bids
from argparse import ArgumentParser
import pandas as pd
import json

from utils import DCM2NII, ADNI_DB, create_dataset_description, init_participants, init_sessions, final_report
from ADNI.setup import *

#################################################################
# This scipt adds new metadata to the image-related JSON files. #
#################################################################

adni_iid_key = 'Image Data ID' #key specifiying adni id on _input_file_
bids_iid_key = 'ImageID' #linked to adni_iid_key

# Specification
datatype = 'anat' # 'pet', 'func'
input_file = '' #csv from ADNI

bids_k = 'ADNIPhase' #column of interest
adni_k = 'COLPROT' #header to add


if not exists(db_file):
    print('Initialising the BIDS dataset.')
    if not exists(bids_dir): makedirs(bids_dir)
    create_dataset_description(bids_dir)

    bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
    bids_db.save(db_file)
else:
    print('Initialising the BIDS dataset. Reading from ' + db_file)
    bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False, database_path=db_file)

input_csv = pd.read_csv(input_file)

subject_list = bids_db.get_subjects()
for sid in subject_list:
    pdb.set_trace()
    json_files = bids_db.get(extension='json', subject=sid, datatype=datatype)
    for jf in json_files:
        json_file = open(jf.path, 'r')
        json_dict = json.load(json_file)
        json_dict[bids_k] = input_csv.loc[input_csv[adni_iid_key] == json_dict[bids_iid_key]].iloc[0][adni_k]
        json_object = json.dumps(json_dict, indent=4)
        json_file = open(jf.path, "w")
        json_file.write(json_object)