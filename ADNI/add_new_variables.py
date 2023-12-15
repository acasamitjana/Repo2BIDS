import pdb
from os.path import join, dirname, exists
from os import makedirs
import bids
from argparse import ArgumentParser
import pandas as pd
import csv
import shutil

from utils import DCM2NII, ADNI_DB, create_dataset_description, init_participants, init_sessions, final_report
from ADNI.setup import *

#############################################################################
# This scipt adds new variables to participants.tsv or *_sessions.tsv files #
#############################################################################

ADNI_PID = 'PTID' #key specifiying adni id on _input_file_
BIDS_PID = 'adni_id' # linked to adni_pid_key

ADNI_SID = 'VISCODE' #key specifiying adni session id on _input_file_
BIDS_SID = 'viscode' # linked to adni_sid_key

# Specification
output_file = 'session' # 'session' or 'participants'
input_file = '/home/acasamitjana/data/util-files/ADNIMERGE.csv' #csv file from ADNI

bids_k = 'HPVol' #header key in the BIDS tsv file.
adni_k = 'Hippocampus' #column of interest from the ADNI csv file.


if not exists(db_file):
    print('Initialising the BIDS dataset.')
    if not exists(bids_dir): makedirs(bids_dir)
    create_dataset_description(bids_dir)

    bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
    bids_db.save(db_file)
else:
    print('Initialising the BIDS dataset. Reading from ' + db_file)
    if exists(db_file_wc): shutil.rmtree(db_file_wc)
    shutil.move(db_file, db_file_wc)
    bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False, database_path=db_file_wc)

input_csv = pd.read_csv(input_file)

p_df = pd.read_csv(join(bids_db.root, 'participants.tsv'), delimiter='\t')
p_dict = {p['participant_id']: p.to_dict() for _, p in p_df.iterrows()}
pdb.set_trace()
if output_file == 'participants':
    to_write = []

    col_names = p_df.columns.to_list() + [bids_k]
    for pid, p_info in p_dict.items():
        # Get rows corresponding to the adni subject
        adni_sbj = input_csv.loc[input_csv[ADNI_PID] == p_info[BIDS_PID]]
        if len(adni_sbj):
            p_info[bids_k] = adni_sbj[adni_k].values[0]

        to_write.append(p_info)

    with open(join(bids_db.root, 'participants.tsv'), 'w') as csv_file:
        tsv_writer = csv.DictWriter(csv_file, fieldnames=col_names, delimiter='\t')
        tsv_writer.writeheader()
        tsv_writer.writerows(to_write)


elif output_file == 'session':

    subject_list = bids_db.get_subjects()
    for pid in subject_list:
        to_write = []
        s_df = pd.read_csv(join(bids_db.root, 'sub-' + pid, 'sub-' + pid + '_sessions.tsv'), delimiter='\t')
        s_dict = {p['session_id']: p.to_dict() for _, p in s_df.iterrows()}
        col_names = s_df.columns.to_list() + [bids_k]
        for sid, s_info in s_dict.items():
            # Get rows corresponding to the adni subject
            adni_sbj = input_csv.loc[input_csv[ADNI_PID] == p_dict[pid][BIDS_PID]]

            # Get rows corresponding to the session
            adni_sess = adni_sbj.loc[input_csv[ADNI_SID] == s_info[BIDS_SID]]
            if len(adni_sess):
                s_info[bids_k] = adni_sess[adni_k].values[0]

            to_write.append(s_info)

        with open(join(bids_db.root, 'sub-' + pid, 'sub-' + pid + '_sessions.tsv'), 'w') as csv_file:
            tsv_writer = csv.DictWriter(csv_file, fieldnames=col_names, delimiter='\t')
            tsv_writer.writeheader()
            tsv_writer.writerows(to_write)

else:
    raise ValueError("Please, provide a valid value for _output_file_: *participants* or *session*.")


if exists(db_file): shutil.rmtree(db_file)
bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
bids_db.save(db_file)
