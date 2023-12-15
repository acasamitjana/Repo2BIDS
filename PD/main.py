############
# Currently, I'm using the NIFTI images directly and filling in the json manually
import pdb
import shutil
from os.path import join, exists, dirname
from os import makedirs, listdir, chmod
import csv
import bids
from datetime import datetime
import numpy as np
import json
import nibabel as nib

IN_ROOT_PATH = '/mnt/nas01/parkinson/'
DOC_PATH = '/home/acasamitjana/data/PD-HCB/Documents/csv'#'/vol/home/student_pd12/ADRIA/Documents/PD_csvs'
BIDS_DIR = '/media/acasamitjana/SEAGATE/PD/PD-BIDS/rawdata'#/vol/home/student_pd12/ADRIA/PD-BIDS/rawdata'
FILE_TREE = {
    '2010': {'t1': {'dir': 'Imatges_originals/MCI_2010/t1_MCI_2010', 'suffix': ''},
             't2': {'dir': 'Imatges_originals/MCI_2010/t2', 'suffix': '_1'},
             'flair': {'dir': 'Imatges_originals/MCI_2010/flair_MCI_2010', 'suffix': '_1'},
             'rest': {'dir': 'Imatges_originals/MCI_2010/rest_MCI_2010', 'suffix': ''}},

    '2013': {'t1': {'dir': 'Imatges_originals/MCI_2013/T1_MCI_2013', 'suffix': '_2'},
             't2': {'dir': 'Imatges_originals/MCI_2013/t2_MCI_2013', 'suffix': '_2'},
             'flair': {'dir': 'Imatges_originals/MCI_2013/FLAIR_MCI_2013', 'suffix': '_2'},
             'rest': {'dir': 'Imatges_originals/MCI_2013/rest_MCI_2013', 'suffix': '_2'}},

    '2015': {'t1': {'dir': 'Imatges_preprocessades/MARATO/2015/reconstruccio_2015/t1', 'suffix': '_3'},
             't2': {'dir': 'Imatges_preprocessades/MARATO/2015/reconstruccio_2015/t2', 'suffix': '_3'},
             'flair': {'dir': 'Imatges_preprocessades/MARATO/2015/reconstruccio_2015/flair', 'suffix': '_3'},
             'rest': {'dir': 'Imatges_preprocessades/MARATO/2015/reconstruccio_2015/rest', 'suffix': '_3'}},

    '2016': {'t1': {'dir': 'Imatges_preprocessades/MARATO/2016/reconstruccio_2016/t1', 'suffix': '_4'},
             't2': {'dir': 'Imatges_preprocessades/MARATO/2016/reconstruccio_2016/t2', 'suffix': '_4'},
             'flair': {'dir': 'Imatges_preprocessades/MARATO/2016/reconstruccio_2016/flair', 'suffix': '_4'},
             'rest': {'dir': 'Imatges_preprocessades/MARATO/2016/reconstruccio_2016/rest', 'suffix': '_4'}},
}

# PATH_2015 = join(IN_ROOT_PATH, 'MCI_2015_Marato_2015_2016', 'LONG_2015_DICOMS', 'LONGITUDINAL_3')
# PATH_2016 = join(IN_ROOT_PATH, 'MCI_2015_Marato_2015_2016', 'LONG_2016_DICOMS')

# Read CSVs de 210, 2013, 2015, 2016
def read_csv(path):
    out_dict = {}
    with open(path, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            out_dict[row['ID']] = row

    return out_dict

# Get a unique list of participants in part_dic = {p1: {}, ...}. For each participant,
# get information related to each session in part_dic[p1] = {sess1: { *metadata* }}

csv_dict = {}
for file in listdir(DOC_PATH):
    csv_dict[file.split('_')[0]] = read_csv(join(DOC_PATH, file))

part_dict = {}
for csv_key, csv_file in csv_dict.items():
    for subj_id, subj_md in csv_file.items():
        if subj_id not in part_dict.keys():
            part_dict[subj_id] = {csv_key: subj_md}
        else:
            part_dict[subj_id][csv_key] = subj_md

pid_list = list(part_dict.keys())

# Insert the associated imaging paths (t1, t2, flair and rs-fMRI). Here I also need to check if it exists
for pid, sess_dict in part_dict.items():
    for sess_id in sess_dict.keys():
        part_dict[pid][sess_id]['t1'] = []
        part_dict[pid][sess_id]['t2'] = []
        part_dict[pid][sess_id]['flair'] = []
        part_dict[pid][sess_id]['rest'] = []


        if exists(join(IN_ROOT_PATH, FILE_TREE[sess_id]['t1']['dir'], pid + FILE_TREE[sess_id]['t1']['suffix'] + '.nii.gz')):
            part_dict[pid][sess_id]['t1'] = [join(IN_ROOT_PATH, FILE_TREE[sess_id]['t1']['dir'], pid + FILE_TREE[sess_id]['t1']['suffix'] + '.nii.gz')]
            if exists(join(IN_ROOT_PATH, FILE_TREE[sess_id]['t1']['dir'], pid + FILE_TREE[sess_id]['t1']['suffix'] + '_2.nii.gz')):
                part_dict[pid][sess_id]['t1'] += [join(IN_ROOT_PATH, FILE_TREE[sess_id]['t1']['dir'], pid + FILE_TREE[sess_id]['t1']['suffix'] + '_2.nii.gz')]

        if exists(join(IN_ROOT_PATH, FILE_TREE[sess_id]['t2']['dir'], pid + FILE_TREE[sess_id]['t2']['suffix'] + '.nii.gz')):
            part_dict[pid][sess_id]['t2'] = [join(IN_ROOT_PATH, FILE_TREE[sess_id]['t2']['dir'], pid + FILE_TREE[sess_id]['t2']['suffix'] + '.nii.gz')]
            if exists(join(IN_ROOT_PATH, FILE_TREE[sess_id]['t2']['dir'], pid + FILE_TREE[sess_id]['t2']['suffix'] + '_2.nii.gz')):
                part_dict[pid][sess_id]['t2'] += [join(IN_ROOT_PATH, FILE_TREE[sess_id]['t2']['dir'], pid + FILE_TREE[sess_id]['t2']['suffix'] + '_2.nii.gz')]

        if exists(join(IN_ROOT_PATH, FILE_TREE[sess_id]['flair']['dir'], pid + FILE_TREE[sess_id]['flair']['suffix'] + '.nii.gz')):
            part_dict[pid][sess_id]['flair'] = [join(IN_ROOT_PATH, FILE_TREE[sess_id]['flair']['dir'], pid + FILE_TREE[sess_id]['flair']['suffix'] + '.nii.gz')]
            if exists(join(IN_ROOT_PATH, FILE_TREE[sess_id]['flair']['dir'], pid + FILE_TREE[sess_id]['flair']['suffix'] + '_2.nii.gz')):
                part_dict[pid][sess_id]['flair'] += [join(IN_ROOT_PATH, FILE_TREE[sess_id]['flair']['dir'], pid + FILE_TREE[sess_id]['flair']['suffix'] + '_2.nii.gz')]

        if exists(join(IN_ROOT_PATH, FILE_TREE[sess_id]['rest']['dir'], pid + FILE_TREE[sess_id]['rest']['suffix'] + '.nii.gz')):
            part_dict[pid][sess_id]['rest'] = [join(IN_ROOT_PATH, FILE_TREE[sess_id]['rest']['dir'], pid + FILE_TREE[sess_id]['rest']['suffix'] + '.nii.gz')]
            if exists(join(IN_ROOT_PATH, FILE_TREE[sess_id]['rest']['dir'], pid + FILE_TREE[sess_id]['rest']['suffix'] + '_2.nii.gz')):
                part_dict[pid][sess_id]['rest'] += [join(IN_ROOT_PATH, FILE_TREE[sess_id]['rest']['dir'], pid + FILE_TREE[sess_id]['rest']['suffix'] + '_2.nii.gz')]

bids_loader = bids.layout.BIDSLayout(root=BIDS_DIR, validate=False)


#
# # Create dataset description
# data_descr = {}
# data_descr['Name'] = 'Parkinson data from Hospital Clinic de Barcelona'
# data_descr['BIDSVersion'] = ''
# data_descr['License'] = ''
# data_descr['Authors'] = ['Barbara Segura']
# data_descr['Acknowledgements'] = ''
# data_descr['HowToAcknowledge'] = ''
# data_descr['Funding'] = ['']
# data_descr['ReferencesAndLinks'] = ['']
# data_descr['DatasetDOI'] = ''
# data_descr_path = join(BIDS_DIR, 'dataset_description.json')
# json_object = json.dumps(data_descr, indent=4)
# with open(data_descr_path, 'w') as outfile:
#     outfile.write(json_object)
#


# Create participants.tsv
p_csv_fields = ['participant_id', 'num_sessions', 'birth_date', 'sex', 'education', 'diagnosis', 'diagnosis_date']

p_csv_dict = []
for pid, sess_dict in part_dict.items():
    for sess_id, sess_md in sess_dict.items():
        p_csv_dict += [{
            'participant_id': pid,
            'num_sessions': len(sess_dict),
            'sex': sess_md['SEX_' + sess_id],
            'birth_date': sess_md['BDAY_' + sess_id],
            'education': sess_md['EDUCATION_yrs_' + sess_id],
            'diagnosis': 'NC' if 'nc' in pid else 'PD',
            'diagnosis_date': sess_md['DIAGNOSIS_DATE_' + sess_id]
        }]
        break

# with open(join(BIDS_DIR, 'participants.tsv'), 'w') as csvfile:
#     csvwriter = csv.DictWriter(csvfile, fieldnames=p_csv_fields, delimiter='\t')
#     csvwriter.writeheader()
#     csvwriter.writerows(p_csv_dict)
#
# chmod(join(BIDS_DIR, 'participants.tsv'), 0o755)

participants_json = {
    "participant_id": {
        "Description": "Participant identifier in BIDS format.",
        },
    "num_sessions": {
        "Description": "Number of sessions per subject at the time of creation.",
        },
    "sex": {
        "Description": "Sex of the participant as reported by the participant",
        "Levels": {
            "0": "male",
            "1": "female"
                }
        },
    "birth_date": {
        "Description": 'Date of birth of the participant',
        },
    "education": {
        "Description": "Education years of the participant",
        "Units": "years"
        },
    "diagnosis ": {
        "Description": "Reported diagnostic category",
        "Levels": {
            "NC": "normal controls",
            "PD": "parkinson's disease participants"
            }
        },
     "diagnosis_date ": {
         "Description": "For Parkinson's disease subjects, year of diagnosis.",
     }
}
# data_descr_path = join(BIDS_DIR, 'participants.json')
# json_object = json.dumps(participants_json, indent=4)
# with open(data_descr_path, 'w') as outfile:
#     outfile.write(json_object)
#
# chmod(data_descr_path, 0o755)

# Create sessions.tsv, COPY images using pybids to build path and fill JSON files.
sess_csv_fields = ['session_id', 'age', 'time_to_bl_days', 'scan_date', 'updrsiii', 'hy', 'disease_duration', 'ledd', 'mmse', 'upsit']
for pid, sess_dict in part_dict.items():
    print('Subject: ' + pid)
    subject_id = 'sub-' + pid
    if not exists(join(BIDS_DIR, subject_id)):
        makedirs(join(BIDS_DIR, subject_id))
        chmod(join(BIDS_DIR, subject_id), 0o755)

    sess_csv_dict = []
    age_days = []
    for sess_y, sess_md in sess_dict.items():
        session_id = 'ses-' + sess_y
        if not exists(join(BIDS_DIR, subject_id, session_id)):
            makedirs(join(BIDS_DIR, subject_id, session_id))
            chmod(join(BIDS_DIR, subject_id, session_id), 0o755)

        # Sessions file
        date_format = "%d/%m/%Y"
        scan_date = datetime.strptime(sess_md['MRI_DATE_' + sess_y], date_format)
        birth_date = datetime.strptime(sess_md['BDAY_' + sess_y], date_format)
        age_days += [(scan_date - birth_date).days]
        age = (scan_date - birth_date).days / 365.25

        sess_csv_dict += [{
            'session_id': sess_y,
            'age': age,
            'scan_date': sess_md['MRI_DATE_' + sess_y],
            'updrsiii': sess_md['UPDRSIII_' + sess_y],
            'hy': sess_md['HY_' + sess_y],
            'disease_duration': sess_md['DISEASE_DURATION_yrs_' + sess_y],
            'ledd': sess_md['LEDD_' + sess_y],
            'mmse': sess_md['MMSE_' + sess_y] if 'MMSE_' + sess_y in sess_md.keys() else '-4',
            'upsit': sess_md['UPSIT_' + sess_y]
        }]

        # # Copy images
        # ## T1
        # for it_im, impath in enumerate(sess_md['t1']):
        #     run_dict = {'run': str(it_im+1).zfill(2)} if len(sess_md['t1']) > 1 else {}
        #     outpath = bids_loader.build_path({**{'subject': pid, 'session': sess_y, 'suffix': 'T1w', 'extension': 'nii.gz'}, **run_dict})
        #     json_path = bids_loader.build_path({**{'subject': pid, 'session': sess_y, 'suffix': 'T1w', 'extension': 'json'}, **run_dict})
        #     if not exists(dirname(outpath)):
        #         makedirs(dirname(outpath))
        #         chmod(dirname(outpath), 0o755)
        #
        #     shutil.copy(impath, outpath)
        #
        #     proxy = nib.load(outpath)
        #     aff = proxy.affine
        #     pixdim = str(np.sqrt(np.sum(aff * aff, axis=0))[:-1])
        #     im_json = {"ScanDate": sess_md['MRI_DATE_' + sess_y],
        #                "Age": str(age),
        #                "VoxelSize": {
        #                    "R": str(pixdim[0]),
        #                    "A": str(pixdim[1]),
        #                    "S": str(pixdim[2])
        #                },
        #                "ImageShape": {
        #                    "X": str(proxy.shape[0]),
        #                    "Y": str(proxy.shape[1]),
        #                    "Z": str(proxy.shape[2])
        #                }}
        #     json_object = json.dumps(im_json, indent=4)
        #     with open(json_path, 'w') as outfile:
        #         outfile.write(json_object)
        #
        #     chmod(outpath, 0o755)
        #     chmod(json_path, 0o755)
        #
        # ## T2
        # for it_im, impath in enumerate(sess_md['t2']):
        #     run_dict = {'run': str(it_im + 1).zfill(2)} if len(sess_md['t2']) > 1 else {}
        #     outpath = bids_loader.build_path( {**{'subject': pid, 'session': sess_y, 'suffix': 'T2w', 'extension': 'nii.gz'}, **run_dict})
        #     json_path = bids_loader.build_path( {**{'subject': pid, 'session': sess_y, 'suffix': 'T2w', 'extension': 'json'}, **run_dict})
        #     if not exists(dirname(outpath)):
        #         makedirs(dirname(outpath))
        #         chmod(dirname(outpath), 0o755)
        #
        #     shutil.copy(impath, outpath)
        #
        #     proxy = nib.load(outpath)
        #     aff = proxy.affine
        #     pixdim = str(np.sqrt(np.sum(aff * aff, axis=0))[:-1])
        #     im_json = {"ScanDate": sess_md['MRI_DATE_' + sess_y],
        #                "Age": str(age),
        #                "VoxelSize": {
        #                    "R": str(pixdim[0]),
        #                    "A": str(pixdim[1]),
        #                    "S": str(pixdim[2])
        #                },
        #                "ImageShape": {
        #                    "X": str(proxy.shape[0]),
        #                    "Y": str(proxy.shape[1]),
        #                    "Z": str(proxy.shape[2])
        #                }}
        #     json_object = json.dumps(im_json, indent=4)
        #     with open(json_path, 'w') as outfile:
        #         outfile.write(json_object)
        #
        #     chmod(outpath, 0o755)
        #     chmod(json_path, 0o755)
        #
        # ## FLAIR
        # for it_im, impath in enumerate(sess_md['flair']):
        #     run_dict = {'run': str(it_im + 1).zfill(2)} if len(sess_md['flair']) > 1 else {}
        #     outpath = bids_loader.build_path( {**{'subject': pid, 'session': sess_y, 'suffix': 'FLAIR', 'extension': 'nii.gz'}, **run_dict})
        #     json_path = bids_loader.build_path( {**{'subject': pid, 'session': sess_y, 'suffix': 'FLAIR', 'extension': 'json'}, **run_dict})
        #     if not exists(dirname(outpath)):
        #         makedirs(dirname(outpath))
        #         chmod(dirname(outpath), 0o755)
        #
        #     shutil.copy(impath, outpath)
        #
        #     proxy = nib.load(outpath)
        #     aff = proxy.affine
        #     pixdim = str(np.sqrt(np.sum(aff * aff, axis=0))[:-1])
        #     im_json = {"ScanDate": sess_md['MRI_DATE_' + sess_y],
        #                "Age": str(age),
        #                "VoxelSize": {
        #                    "R": str(pixdim[0]),
        #                    "A": str(pixdim[1]),
        #                    "S": str(pixdim[2])
        #                },
        #                "ImageShape": {
        #                    "X": str(proxy.shape[0]),
        #                    "Y": str(proxy.shape[1]),
        #                    "Z": str(proxy.shape[2])
        #                }}
        #     json_object = json.dumps(im_json, indent=4)
        #     with open(json_path, 'w') as outfile:
        #         outfile.write(json_object)
        #
        #     chmod(outpath, 0o755)
        #     chmod(json_path, 0o755)
        #
        # ## Rest
        # for it_im, impath in enumerate(sess_md['rest']):
        #     run_dict = {'run': str(it_im + 1).zfill(2)} if len(sess_md['rest']) > 1 else {}
        #     outpath = bids_loader.build_path({**{'subject': pid, 'session': sess_y, 'task': 'rest','suffix': 'bold', 'extension': 'nii.gz'}, **run_dict})
        #     json_path = bids_loader.build_path({**{'subject': pid, 'session': sess_y, 'task': 'rest','suffix': 'bold', 'extension': 'json'}, **run_dict})
        #     if not exists(dirname(outpath)):
        #         makedirs(dirname(outpath))
        #         chmod(dirname(outpath), 0o755)
        #     shutil.copy(impath, outpath)
        #
        #     proxy = nib.load(outpath)
        #     aff = proxy.affine
        #     pixdim = str(np.sqrt(np.sum(aff * aff, axis=0))[:-1])
        #     im_json = {"ScanDate": sess_md['MRI_DATE_' + sess_y],
        #                "Age": str(age),
        #                "VoxelSize": {
        #                    "R": str(pixdim[0]),
        #                    "A": str(pixdim[1]),
        #                    "S": str(pixdim[2])
        #                },
        #                "ImageShape": {
        #                    "X": str(proxy.shape[0]),
        #                    "Y": str(proxy.shape[1]),
        #                    "Z": str(proxy.shape[2])
        #                },
        #                "NumVolumes": str(proxy.shape[3])
        #                }
        #     json_object = json.dumps(im_json, indent=4)
        #     with open(json_path, 'w') as outfile:
        #         outfile.write(json_object)
        #
        #     chmod(outpath, 0o755)
        #     chmod(json_path, 0o755)

    age_bl_days = np.min(age_days)
    for it_sess in range(len(sess_csv_dict)):
        sess_csv_dict[it_sess]['time_to_bl_days'] = age_days[it_sess] - age_bl_days

    with open(join(BIDS_DIR, subject_id, subject_id + '_sessions.tsv'), 'w') as csvfile:
        csvwriter = csv.DictWriter(csvfile, fieldnames=sess_csv_fields, delimiter='\t')
        csvwriter.writeheader()
        csvwriter.writerows(sess_csv_dict)

    chmod(join(BIDS_DIR, subject_id, subject_id + '_sessions.tsv'), 0o755)

