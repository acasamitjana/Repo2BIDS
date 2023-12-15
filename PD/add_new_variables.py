import pdb
from os.path import join, dirname, exists
from os import makedirs
import bids
from argparse import ArgumentParser
import pandas as pd
import csv
import shutil

#############################################################################
# This scipt adds new variables to participants.tsv or *_sessions.tsv files #
#############################################################################

PD_PID = 'ID' #key specifiying adni id on _input_file_
BIDS_PID = 'participant_id' # linked to adni_pid_key

PD_SID = '' #key specifiying adni session id on _input_file_
BIDS_SID = '' # linked to adni_sid_key

# Specification
output_file = 'participants' # 'session' or 'participants'
input_file = '/media/biofisica/BIG_DATA/PD/2010_RM_Cluster_MCI_variables.xlsx' #csv file from ADNI

bids_k = 'Cluster_2010' #header key in the BIDS tsv file.
pd_k = 'Cluster' #column of interest from the ADNI csv file.

# Path
bids_dir = '/media/biofisica/BIG_DATA/PD/PD-BIDS/rawdata'

# if not exists(db_file):
#     print('Initialising the BIDS dataset.')
#     if not exists(bids_dir): makedirs(bids_dir)
#     create_dataset_description(bids_dir)
#
#     bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
#     bids_db.save(db_file)
# else:
#     print('Initialising the BIDS dataset. Reading from ' + db_file)
#     if exists(db_file_wc): shutil.rmtree(db_file_wc)
#     shutil.move(db_file, db_file_wc)
#     bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False, database_path=db_file_wc)

bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
input_csv = pd.read_excel(input_file, sheet_name=0)

p_df = pd.read_csv(join(bids_db.root, 'participants.tsv'), delimiter='\t')
p_dict = {p['participant_id']: p.to_dict() for _, p in p_df.iterrows()}
if output_file == 'participants':
    to_write = []

    col_names = p_df.columns.to_list() + [bids_k]
    for pid, p_info in p_dict.items():
        # Get rows corresponding to the adni subject
        adni_sbj = input_csv.loc[input_csv[PD_PID] == p_info[BIDS_PID]]
        if len(adni_sbj):
            p_info[bids_k] = adni_sbj[pd_k].values[0]

        to_write.append(p_info)

    with open(join(bids_db.root, 'participants.tsv'), 'w') as csv_file:
        tsv_writer = csv.DictWriter(csv_file, fieldnames=col_names, delimiter='\t')
        tsv_writer.writeheader()
        tsv_writer.writerows(to_write)


elif output_file == 'session':

    subject_list = bids_db.get_subjects()
    for pid in subject_list:
        to_write = []
        pdb.set_trace()
        s_df = pd.read_csv(join(bids_db.root, 'sub-' + pid, 'sub-' + pid + '_sessions.tsv'), delimiter='\t')
        s_dict = {p['session_id']: p.to_dict() for _, p in s_df.iterrows()}
        col_names = s_df.columns.to_list() + [bids_k]
        for sid, s_info in s_dict.items():
            # Get rows corresponding to the adni subject
            adni_sbj = input_csv.loc[input_csv[PD_PID] == p_dict[pid][BIDS_PID]]

            # Get rows corresponding to the session
            adni_sess = adni_sbj.loc[input_csv[PD_SID] == s_info[BIDS_SID]]
            if len(adni_sess):
                s_info[bids_k] = adni_sess[pd_k].values[0]

            to_write.append(s_info)

        with open(join(bids_db.root, 'sub-' + pid, 'sub-' + pid + '_sessions.tsv'), 'w') as csv_file:
            tsv_writer = csv.DictWriter(csv_file, fieldnames=col_names, delimiter='\t')
            tsv_writer.writeheader()
            tsv_writer.writerows(to_write)

else:
    raise ValueError("Please, provide a valid value for _output_file_: *participants* or *session*.")


# if exists(db_file): shutil.rmtree(db_file)
# bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
# bids_db.save(db_file)
