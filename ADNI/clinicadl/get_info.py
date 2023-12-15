import csv
import pdb
from os.path import join, dirname, exists
import bids
import nibabel as nib
import numpy as np
from matplotlib import pyplot as plt, font_manager as fm
import time
import pandas as pd

# Input variables
# path_data = '/media/biofisica/BIG_DATA/'
# path_utils = '/home/biofisica/data/util-files'
# bids_dir = '/media/biofisica/BIG_DATA/ADNI-BIDS'

################
# Let's start! #
################
# Read target subjects:
cn_dict = {}
mci_dict = {}
ad_dict = {}
with open('data/CN_train_baseline.tsv', 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter='\t')
    for row in csvreader:
        ses_id = row['session_id'].split('-')[-1][-2:]
        if ses_id not in cn_dict.keys():
            cn_dict[row['participant_id']] = []
        cn_dict[row['participant_id']].append(ses_id)

with open('data/AD_train_baseline.tsv', 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter='\t')
    for row in csvreader:
        ses_id = row['session_id'].split('-')[-1][-2:]
        if ses_id not in ad_dict.keys():
            ad_dict[row['participant_id']] = []
        ad_dict[row['participant_id']].append(ses_id)

with open('data/MCI_baseline.tsv', 'r') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter='\t')
    for row in csvreader:
        ses_id = row['session_id'].split('-')[-1][-2:]
        if ses_id not in mci_dict.keys():
            mci_dict[row['participant_id']] = []
        mci_dict[row['participant_id']].append(ses_id)

data_dict = mci_dict#{**cn_dict, **ad_dict}

info_dict = {}
info_dict['subject_id'] = []
info_dict['session_id'] = []
info_dict['tiv'] = []
info_dict['fld_str'] = []
info_dict['adniprot'] = []
info_dict['age'] = []
info_dict['qualifications'] = []
info_dict['apoe4'] = []
info_dict['sex'] = []
info_dict['ethnicity'] = []
info_dict['race'] = []
info_dict['dx'] = []


subject_dict = {}
adni_merge = pd.read_csv('/home/acasamitjana/Data/ADNI/ADNIMERGE.csv')
adni_merge = adni_merge.set_index(['PTID', 'Month'])
for subject, session_list in data_dict.items():
    sid = '_S_'.join(subject[8:].split('S'))
    for sess in session_list:
        #
        session_data = adni_merge.loc[sid].loc[int(sess)]

        adniprot = session_data['COLPROT']
        try:
            fld_strength = float(session_data['FLDSTRENG'][0])
        except:
            fld_strength = np.nan
            print(session_data['FLDSTRENG'])

        tiv = float(session_data['WholeBrain']) * 1e-6

        age = session_data['AGE']
        qualifications = session_data['PTEDUCAT']
        apoe4 = session_data['APOE4']
        sex = session_data['PTGENDER']
        ethnicity = session_data['PTETHCAT']
        race = session_data['PTRACCAT']
        dx = session_data['DX_bl']

        subject_dict[subject] = {'subject_id': subject, 'session_id': 'ses-M'+sess,  'tiv': tiv, 'fld_str': fld_strength, 'adniprot': adniprot}

        info_dict['subject_id'] += [subject]
        info_dict['session_id'] += ['ses-M'+sess]
        info_dict['tiv'] += [tiv]
        info_dict['fld_str'] += [fld_strength]
        info_dict['adniprot'] += [adniprot]

        info_dict['age'] += [age]
        info_dict['qualifications'] += [qualifications]
        info_dict['apoe4'] += [apoe4]
        info_dict['sex'] += [sex]
        info_dict['ethnicity'] += [ethnicity]
        info_dict['race'] += [race]
        info_dict['dx'] += [dx]

pd.DataFrame.from_dict(info_dict).to_csv('info_mci.csv', sep=' ')

# # Reading BIDS dataset
# db_file =join(dirname(bids_dir), 'ADNI-BIDS_20230616_prep.db')
# if not exists(db_file):
#     print('Initialising the BIDS dataset.')
#     bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
#     bids_db.add_derivatives(join(bids_dir, 'derivatives', 'synthseg'))
#     bids_db.add_derivatives(join(bids_dir, 'derivatives', 'synthsr'))
#     bids_db.save(db_file)
# else:
#     print('Initialising the BIDS dataset. Reading from ' + db_file)
#     bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False, database_path=db_file)
#
# for subject in cn_dict.keys():
#     t1w_image = bids_db.get(subject=subject[4:], session='m000', extension='.nii.gz', suffix='T1w')
#     pdb.set_trace()