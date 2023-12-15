import csv
import pdb
from os import listdir, makedirs
from os.path import join, dirname, exists
import shutil

import nibabel as nib
import numpy as np
from matplotlib import pyplot as plt, font_manager as fm
import time
import pandas as pd

# Input variables
bids_dir = '/mnt/HDD/Data/ADNI/rawdata'
output_dir = '/mnt/HDD/Data/ADNI/clinicadl-validation/rawdata'
output_mci_dir = '/mnt/HDD/Data/ADNI/clinicadl-validation-mci/rawdata'

################
# Let's start! #
################
# Read target subjects:
for it_split in range(5):
    print('Split ' + str(it_split))
    cn_dict = {}
    ad_dict = {}

    with open('data/validation-splits/split-' + str(it_split) +  '/CN_baseline.tsv', 'r') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter='\t')
        for row in csvreader:
            ses_id = row['session_id'].split('-')[-1][-2:]
            cn_dict[row['participant_id']] = {'session': 'ses-m' + ses_id, 'orig_path': None, 'final_path': None}

    with open('data/validation-splits/split-' + str(it_split) + '/AD_baseline.tsv', 'r') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter='\t')
        for row in csvreader:
            ses_id = row['session_id'].split('-')[-1][-2:]
            ad_dict[row['participant_id']] = {'session': 'ses-m' + ses_id, 'orig_path': None, 'final_path': None}

    for pid, pdict in cn_dict.items():
        image_dir = join(bids_dir, pid, pdict['session'], 'anat')
        output_image_dir = join(output_dir, pid, pdict['session'], 'anat')
        if not exists(output_image_dir): makedirs(output_image_dir)

        if not exists(image_dir):
            print(pid)
        else:
            image_path = list(filter(lambda x: 'T1w' in x and 'nii.gz' in x, listdir(image_dir)))
            cn_dict[pid] = {'orig_path': join(image_dir, image_path[0]), 'final_path': join(output_image_dir, image_path[0])}
            shutil.copy(join(image_dir, image_path[0]), join(output_image_dir, image_path[0]))

    for pid, pdict in ad_dict.items():
        image_dir = join(bids_dir, pid, pdict['session'], 'anat')
        output_image_dir = join(output_dir, pid, pdict['session'], 'anat')
        if not exists(output_image_dir): makedirs(output_image_dir)

        if not exists(image_dir):
            print(pid)
        else:
            image_path = list(filter(lambda x: 'T1w' in x and 'nii.gz' in x, listdir(image_dir)))
            ad_dict[pid] = {'orig_path': join(image_dir, image_path[0]), 'final_path': join(output_image_dir, image_path[0])}
            shutil.copy(join(image_dir, image_path[0]), join(output_image_dir, image_path[0]))


for it_split in range(5):
    print('Split ' + str(it_split))
    pmci_dict = {}
    smci_dict = {}
    with open('data/validation-splits/split-' + str(it_split) +  '/pMCI_baseline.tsv', 'r') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter='\t')
        for row in csvreader:
            ses_id = row['session_id'].split('-')[-1][-2:]
            pmci_dict[row['participant_id']] = {'session': 'ses-m' + ses_id, 'orig_path': None, 'final_path': None}

    with open('data/validation-splits/split-' + str(it_split) + '/sMCI_baseline.tsv', 'r') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter='\t')
        for row in csvreader:
            ses_id = row['session_id'].split('-')[-1][-2:]
            smci_dict[row['participant_id']] = {'session': 'ses-m' + ses_id, 'orig_path': None, 'final_path': None}

    for pid, pdict in smci_dict.items():
        image_dir = join(bids_dir, pid, pdict['session'], 'anat')
        output_image_dir = join(output_mci_dir, pid, pdict['session'], 'anat')
        if not exists(output_image_dir): makedirs(output_image_dir)

        if not exists(image_dir):
            print(pid + '_' + pdict['session'])
        else:
            image_path = list(filter(lambda x: 'T1w' in x and 'nii.gz' in x, listdir(image_dir)))
            smci_dict[pid] = {'orig_path': join(image_dir, image_path[0]), 'final_path': join(output_image_dir, image_path[0])}
            shutil.copy(join(image_dir, image_path[0]), join(output_image_dir, image_path[0]))

    for pid, pdict in pmci_dict.items():
        image_dir = join(bids_dir, pid, pdict['session'], 'anat')
        output_image_dir = join(output_mci_dir, pid, pdict['session'], 'anat')
        if not exists(output_image_dir): makedirs(output_image_dir)

        if not exists(image_dir):
            print(pid + '_' + pdict['session'])
        else:
            image_path = list(filter(lambda x: 'T1w' in x and 'nii.gz' in x, listdir(image_dir)))
            pmci_dict[pid] = {'orig_path': join(image_dir, image_path[0]), 'final_path': join(output_image_dir, image_path[0])}
            shutil.copy(join(image_dir, image_path[0]), join(output_image_dir, image_path[0]))
