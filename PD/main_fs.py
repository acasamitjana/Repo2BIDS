############
# Currently, I'm using the NIFTI images directly and filling in the json manually
import pdb
import shutil
from os.path import join, exists, dirname, basename
from os import makedirs, listdir, chmod, environ
import csv
import bids
from datetime import datetime
import numpy as np
import json
import nibabel as nib

DOC_PATH = '/vol/home/student_pd12/ADRIA/Documents/PD_csvs'
BIDS_DIR = '/vol/home/student_pd12/ADRIA/PD-BIDS/rawdata'
FS_DIR = '/vol/home/student_pd12/ADRIA/PD-BIDS/derivatives/freesurfer_2013'
FILE_TREE = {
    # '2010': {'name': '/mnt/nas01/parkinson/Imatges_preprocessades/recon_all_all_projects_NEW_QCACHE/PD_PRE_CTH_2010', 'suffix': ''},
    '2013':  {'name': '/mnt/nas01/parkinson/Imatges_preprocessades/recon_all_all_projects_NEW_QCACHE/PD_POST_CTH_2013', 'suffix': '_2'},
    # '2015':  {'name': '/mnt/nas01/parkinson/Imatges_preprocessades/recon_all_all_projects_NEW_QCACHE/RECON_ALL_Marato_2015', 'suffix': '_3'},
    # '2016':  {'name': '/mnt/nas01/parkinson/Imatges_preprocessades/recon_all_all_projects_NEW_QCACHE/RECON_ALL_Marato_2016', 'suffix': '_4'},
}

# Read CSVs de 210, 2013, 2015, 2016
def read_csv(path):
    out_dict = {}
    with open(path, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            out_dict[row['ID']] = row

    return out_dict

def read_FS_volumes(file):
    etiv = 0
    fs_vols = {}
    start = False
    with open(file, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=' ')
        for row in csvreader:
            row_cool = list(filter(lambda x: x != '', row))
            if start is True:
                fs_vols[int(row_cool[1])] = {'name': row_cool[4], 'vol': float(row_cool[3])}

            if 'ColHeaders' in row_cool and start is False:
                start = True
            elif 'ICV' in row_cool:
                etiv = float(row_cool[-2].split(',')[0])

    # vols = {**{lab: float(fs_vols[lab_str]) for lab, lab_str in labs.items() if 'Thalamus' not in lab_str},
    #         **{lab: float(fs_vols[lab_str + '-Proper']) for lab, lab_str in labs.items() if 'Thalamus' in lab_str}}

    return fs_vols, etiv

fs_lut = {'Background': 0}
with open(join(environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt'), 'r') as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
        info = [r for r in row[None][0].split(' ') if r != '']
        if len(info) < 5: continue
        try:
            fs_lut[info[1]] = int(info[0])
        except:
            continue

def read_FS_aparc(file, fsl_lut):
    etiv = 0
    fs_vols = {}
    start = False
    hemi = basename(file).split('.')[0]
    assert hemi in ['rh', 'lh']

    with open(file, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=' ')
        for row in csvreader:
            row_cool = list(filter(lambda x: x != '', row))
            if start is True:
                lab = fsl_lut['ctx-' + hemi + '-' + row_cool[0]]
                fs_vols[lab] = {'name': row_cool[0], 'area': float(row_cool[2]), 'vol': float(row_cool[3]), 'th': float(row_cool[4]), 'th_std': float(row_cool[5])}

            if 'ColHeaders' in row_cool and start is False:
                start = True

    # vols = {**{lab: float(fs_vols[lab_str]) for lab, lab_str in labs.items() if 'Thalamus' not in lab_str},
    #         **{lab: float(fs_vols[lab_str + '-Proper']) for lab, lab_str in labs.items() if 'Thalamus' in lab_str}}

    return fs_vols




# Get a unique list of participants in part_dic = {p1: {}, ...}. For each participant,
# get information related to each session in part_dic[p1] = {sess1: { *metadata* }}

csv_dict = {}
for file in listdir(DOC_PATH):
    csv_dict[file.split('_')[0]] = read_csv(join(DOC_PATH, file))

part_dict = {}
for csv_key, csv_file in csv_dict.items():
    if csv_key not in FILE_TREE: continue
    for subj_id, subj_md in csv_file.items():
        if subj_id not in part_dict.keys():
            part_dict[subj_id] = {csv_key: subj_md}
        else:
            part_dict[subj_id][csv_key] = subj_md

pid_list = list(part_dict.keys())

# Insert the associated imaging paths (t1, t2, flair and rs-fMRI). Here I also need to check if it exists
missing_keys = {}
for pid, sess_dict in part_dict.items():
    for sess_id in sess_dict.keys():
        pid_fs_dir = join(FILE_TREE[sess_id]['name'], pid + FILE_TREE[sess_id]['suffix'])
        if not exists(pid_fs_dir):

            potential_matches = list(filter(lambda x: pid in x and len(x.split('_')) == 2, listdir(FILE_TREE[sess_id]['name'])))
            if len(potential_matches) > 0:
                pid_fs_dir = join(FILE_TREE[sess_id]['name'], potential_matches[0])
            else:
                if pid not in missing_keys.keys(): missing_keys[pid] = []
                missing_keys[pid] += [sess_id]
                continue
        part_dict[pid][sess_id]['t1'] = join(pid_fs_dir, 'mri', 'nu.mgz')
        part_dict[pid][sess_id]['seg'] = join(pid_fs_dir, 'mri', 'aparc+aseg.mgz')

        aseg_stats = join(pid_fs_dir, 'stats', 'aseg.stats')
        lh_aparc_stats = join(pid_fs_dir, 'stats', 'lh.aparc.stats')
        rh_aparc_stats = join(pid_fs_dir, 'stats', 'rh.aparc.stats')
        aseg_vols, etiv = read_FS_volumes(aseg_stats)
        lh_aparc_vols = read_FS_aparc(lh_aparc_stats, fs_lut)
        rh_aparc_vols = read_FS_aparc(rh_aparc_stats, fs_lut)

        part_dict[pid][sess_id]['etiv'] = etiv
        part_dict[pid][sess_id]['volume_dict'] = {
            **{k: v['vol'] for k, v in aseg_vols.items()},
            **{k: v['vol'] for k, v in lh_aparc_vols.items()},
            **{k: v['vol'] for k, v in rh_aparc_vols.items()}
        }

        part_dict[pid][sess_id]['thickness_dict'] = {
            **{k: v['th'] for k, v in lh_aparc_vols.items()},
            **{k: v['th'] for k, v in rh_aparc_vols.items()}
        }

        part_dict[pid][sess_id]['area_dict'] = {
            **{k: v['area'] for k, v in lh_aparc_vols.items()},
            **{k: v['area'] for k, v in rh_aparc_vols.items()}
        }

for pid, sess_list in missing_keys.items():
    for s in sess_list:
        part_dict[pid].pop(s)

bids_loader = bids.layout.BIDSLayout(root=BIDS_DIR, validate=False)

# Create sessions.tsv, COPY images using pybids to build path and fill JSON files.
for pid, sess_dict in part_dict.items():
    print('Subject: ' + pid)
    subject_id = 'sub-' + pid
    if not exists(join(FS_DIR, subject_id)):
        makedirs(join(FS_DIR, subject_id))
        chmod(join(FS_DIR, subject_id), 0o755)


    for sess_y, sess_md in sess_dict.items():
        session_id = 'ses-' + sess_y
        if not exists(join(FS_DIR, subject_id, session_id)):
            makedirs(join(FS_DIR, subject_id, session_id))
            chmod(join(FS_DIR, subject_id, session_id), 0o755)

        # Copy images
        ## T1
        impath = basename(bids_loader.build_path({'subject': pid, 'session': sess_y, 'suffix': 'T1w', 'extension': 'nii.gz'}))
        imjson_path = impath.replace('nii.gz', 'json')
        segpath = impath.replace('T1w', 'dseg')
        segjson_path = segpath.replace('nii.gz', 'json')
        segvol_path = segjson_path.replace('json', 'tsv')
        segarea_path = segvol_path.replace('dseg', 'area')
        segth_path = segvol_path.replace('dseg', 'thickness')

        img = nib.load(sess_md['t1'])
        nib.save(img, join(FS_DIR, subject_id, session_id, impath))
        aff = img.affine
        pixdim = str(np.sqrt(np.sum(aff * aff, axis=0))[:-1])
        im_json = {"ScanDate": sess_md['MRI_DATE_' + sess_y],
                   "VoxelSize": {
                       "R": str(pixdim[0]),
                       "A": str(pixdim[1]),
                       "S": str(pixdim[2])
                   },
                   "ImageShape": {
                       "X": str(img.shape[0]),
                       "Y": str(img.shape[1]),
                       "Z": str(img.shape[2])
                   }}
        json_object = json.dumps(im_json, indent=4)
        with open(join(FS_DIR, subject_id, session_id, imjson_path), 'w') as outfile:
            outfile.write(json_object)

        img = nib.load(sess_md['seg'])
        nib.save(img, join(FS_DIR, subject_id, session_id, segpath))
        aff = img.affine
        pixdim = str(np.sqrt(np.sum(aff * aff, axis=0))[:-1])
        im_json = {"VoxelSize": {
                       "R": str(pixdim[0]),
                       "A": str(pixdim[1]),
                       "S": str(pixdim[2])
                   },
                   "ImageShape": {
                       "X": str(img.shape[0]),
                       "Y": str(img.shape[1]),
                       "Z": str(img.shape[2])
                   }}
        json_object = json.dumps(im_json, indent=4)
        with open(join(FS_DIR, subject_id, session_id, segjson_path), 'w') as outfile:
            outfile.write(json_object)


        with open(join(FS_DIR, subject_id, session_id, segvol_path), 'w') as csvfile:
            csvwriter = csv.DictWriter(csvfile, fieldnames= list(sess_md['volume_dict'].keys()), delimiter='\t')
            csvwriter.writeheader()
            csvwriter.writerow(sess_md['volume_dict'])

        with open(join(FS_DIR, subject_id, session_id, segth_path), 'w') as csvfile:
            csvwriter = csv.DictWriter(csvfile, fieldnames= list(sess_md['thickness_dict'].keys()), delimiter='\t')
            csvwriter.writeheader()
            csvwriter.writerow(sess_md['thickness_dict'])

        with open(join(FS_DIR, subject_id, session_id, segarea_path), 'w') as csvfile:
            csvwriter = csv.DictWriter(csvfile, fieldnames= list(sess_md['area_dict'].keys()), delimiter='\t')
            csvwriter.writeheader()
            csvwriter.writerow(sess_md['area_dict'])

        np.save(join(FS_DIR, subject_id, session_id, 'etiv.npy'), sess_md['etiv'])
        chmod(join(FS_DIR, subject_id, session_id, segvol_path), 0o755)
        chmod(join(FS_DIR, subject_id, session_id, segth_path), 0o755)
        chmod(join(FS_DIR, subject_id, session_id, segarea_path), 0o755)
        chmod(join(FS_DIR, subject_id, session_id, 'etiv.npy'), 0o755)


if not exists(join(FS_DIR, 'freesurfer_lut.txt')):
    ASEG_DICT = {
        0: 'background',
        2: 'left cerebral white matter',
        3: 'left cerebral cortex',
        4: 'left lateral ventricle',
        5: 'left inferior lateral ventricle',
        7: 'left cerebellum white matter',
        8: 'left cerebellum cortex',
        10: 'left thalamus',
        11: 'left caudate',
        12: 'left putamen',
        13: 'left pallidum',
        14: '3rd ventricle',
        15: '4th ventricle',
        16: 'brain-stem',
        17: 'left hippocampus',
        18: 'left amygdala',
        24: 'CSF',
        26: 'left accumbens area',
        28: 'left ventral DC',
        30: 'left vessel',
        31: 'left choroid-plexus',
        41: 'right cerebral white matter',
        42: 'right cerebral cortex',
        43: 'right lateral ventricle',
        44: 'right inferior lateral ventricle',
        46: 'right cerebellum white matter',
        47: 'right cerebellum cortex',
        49: 'right thalamus',
        50: 'right caudate',
        51: 'right putamen',
        52: 'right pallidum',
        53: 'right hippocampus',
        54: 'right amygdala',
        58: 'right accumbens area',
        60: 'right ventral DC',
        62: 'right vessel',
        63: 'right choroid plexus',
        77: 'WM hypo',
        80: 'non WM hypo',
        85: 'optic chiasm',
        251: 'cc posterior',
        252: 'cc mid posterior',
        253: 'cc central',
        254: 'cc mid anterior',
        255: 'cc anterior'
    }
    ASEG_DICT = {k: ASEG_DICT[k] for k in sorted(ASEG_DICT.keys(), key=lambda x: x)}
    ASEG_DICT_REV = {v: k for k, v in ASEG_DICT.items()}
    ASEG_LUT = {k: it_k for it_k, k in enumerate(ASEG_DICT.keys())}
    ASEG_ARR = np.array(list(ASEG_DICT.keys()), dtype='int')
    APARC_ARR = np.concatenate((np.arange(1, 4) + 1000, np.arange(5, 36) + 1000, np.arange(1, 4) + 2000, np.arange(5, 36) + 2000), axis=0)
    ASEG_APARC_ARR = np.unique(np.concatenate((ASEG_ARR, APARC_ARR), axis=0))

    # --------------- #

    labels_abbr = {
        0: 'BG',
        2: 'L-Cerebral-WM',
        3: 'L-Cerebral-GM',
        4: 'L-Lat-Vent',
        5: 'L-Inf-Lat-Vent',
        7: 'L-Cerebell-WM',
        8: 'L-Cerebell-GM',
        10: 'L-TH',
        11: 'L-CAU',
        12: 'L-PU',
        13: 'L-PA',
        14: '3-Vent',
        15: '4-Vent',
        16: 'BS',
        17: 'L-HIPP',
        18: 'L-AM',
        26: 'L-ACC',
        28: 'L-VDC',
        41: 'R-Cerebral-WM',
        42: 'R-Cerebral-GM',
        43: 'R-Lat-Vent',
        44: 'R-Inf-Lat-Vent',
        46: 'R-Cerebell-WM',
        47: 'R-Cerebell-WM',
        49: 'R-TH',
        50: 'R-CAU',
        51: 'R-PU',
        52: 'R-PA',
        53: 'R-HIPP',
        54: 'R-AM',
        58: 'R-ACC',
        60: 'R-VDC',
    }

    fs_lut = {0: {'name': 'Background', 'R': 0, 'G': 0, 'B': 0}}
    with open(join(environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt'), 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            info = [r for r in row[None][0].split(' ') if r != '']
            if len(info) < 5: continue
            try:
                fs_lut[int(info[0])] = {'name': info[1], 'R': info[2], 'G': info[3], 'B': info[4]}
            except:
                continue

    header = ['index', 'name', 'abbreviation', 'R', 'G', 'B', 'mapping']
    label_dict = [
        {'index': label, 'name': fs_lut[label]['name'],
         'abbreviation': labels_abbr[label] if label in labels_abbr else fs_lut[label]['name'],
         'R': fs_lut[label]['R'], 'G': fs_lut[label]['G'], 'B': fs_lut[label]['B'], 'mapping': it_label}
        for it_label, label in enumerate(ASEG_APARC_ARR)
    ]

    with open(join(FS_DIR, 'freesurfer_lut.txt'), 'w') as csvfile:
        csvreader = csv.DictWriter(csvfile, fieldnames=header, delimiter='\t')
        csvreader.writeheader()
        csvreader.writerows(label_dict)

print('\nDone\n')

