import pdb
from os.path import join, isdir, dirname, basename, exists
from os import listdir, makedirs, rename
from datetime import datetime
import subprocess
import json
from joblib import delayed, Parallel
import sys
import csv
import shutil

from zipfile import ZipFile
import numpy as np
import pandas as pd
from ADNI.setup import *


def replace_underscore(word):
    return word.replace('_', '')

def final_report(sample_df, sample_filter_df, failed_images):
    with open('final_output.txt', 'w') as txtfile:
        sys.stdout = txtfile  # Change the standard output to the file we created.
        print('\n\n\n####################')
        print('##### Summary ######')
        print('####################')
        print('\n Total Images: ' + str(len(sample_df)))
        print('\n Total Filtered Images: ' + str(len(sample_filter_df)))
        nifti_paths = {image_df['dicom_path']: image_df['nifti_path'] for _, image_df in sample_filter_df.iterrows()}
        rep_im_d = {}
        for d, n in nifti_paths.items():
            if n in rep_im_d.keys():
                rep_im_d[n].append(d)
            else:
                rep_im_d[n] = [d]
        rep_im_nifti = [n for n, dlist in rep_im_d.items() if len(dlist) > 1]

        print('\n -'.join([image_df['dicom_path'] for _, image_df in sample_filter_df.iterrows()]))
        print('\nTotal Converted Images: ' + str(len(sample_filter_df) - len(failed_images)))
        print('\nTotal Failed Images: ' + str(len(failed_images)))
        print('\n -'.join([fi['dicom_path'] + ': ' + fi['error'] for fi in failed_images]))
        print('\nTotal Repeated Images: ' + str(len(rep_im_nifti)))
        print('\n -' + '\n -'.join([n + ': ' + '\n   * ' + '\n   * '.join(rep_im_d[n]) for n in rep_im_nifti]))


class ADNI_DB(object):
    '''
    The structure of the downloaded dataset is the following:
    * path_data: directory containing the downloaded data. At least, there should be two directories:
        - ADNI-zips: all zip files from ADNI containing the images.
        - ADNI-csv: the csv files related to each zip file downloaded. They contain information like Modality and
                    Description. This CSV is automatically downloaded from ADNI with each download.
    * path_utils: directory containing Study Data downloaded from ADNI. It should contain, at least, the MRILIST,
                  the PET_META_LIST, the VISITS and ADNIMERGE files.

    The output is two directories:
    * path_data/ADNI-wc: the working temporal directory where one should unzip the required ADNI files.
    * path_data/ADNI-BIDS: the dataset transformed to BIDS.

    '''
    def __init__(self, path_data, path_utils, bids_db, load_df=True):
        self.path_data = path_data
        self.path_utils = path_utils
        self.load_df = load_df
        self.bids_db = bids_db

        # Extra fields added to the columns of self.sample_df_filepath (i.e., downloaded CSVs associated to each ZIP)
        self.extra_download_info = ['dicom_path', 'nifti_path', 'subject_id', 'session_id', 'file_dict']


    def build(self):
        # Working directory to unzip all files
        self.working_dir = join(self.path_data, 'ADNI-wc')

        # Directory containing zip files downloaded
        self.zips_dir = join(self.path_data, 'ADNI-zips')
        # Directory containing csv files downloaded linked to each zip file
        self.csvs_dir = join(self.path_data, 'ADNI-csvs')
        # CSV file as a result of merging all downloaded csv. Used to read all downloaded files before unzipping them.
        self.sample_df_filepath = join(self.path_utils, 'sample_df.csv')

        # MRILIST file merged with VISITS File
        self.mrilist_filepath = join(self.path_utils, 'MRILIST-with-Visit.csv')
        # PET_META_LIST file merged with VISITS File
        self.petlist_filepath = join(self.path_utils, 'PET_META_LIST-with-Visit.csv')
        # VISITS file
        self.visits_filepath = join(self.path_utils, 'VISITS.csv')
        # ADNIMERGE file
        self.adnimerge_filepath = join(self.path_utils, 'ADNIMERGE.csv')

        # All zips files
        self.zips_download = listdir(self.zips_dir)

        self.mri_df = pd.read_csv(self.mrilist_filepath)
        self.pet_df = pd.read_csv(self.petlist_filepath)
        self.visits_df = pd.read_csv(self.visits_filepath)
        self.adnimerge_df = pd.read_csv(self.adnimerge_filepath, low_memory=False)

    def read_csvs(self, csv_files):
        '''
        This function is used if one wants to read a specific CSV file (instead of all files found in self.csvs_dir).
        It could be useful if a new extra package is downloaded.

        :param csv_files: list of files (filenames which should be stored in self.csvs_dir).
        :return: pd.Dataframe containing a single table for all downloads.
        '''
        print(' o Reading csv files.')

        files = []
        for filename in csv_files:
            if filename[0] == '.': continue
            df = pd.read_csv(join(self.csvs_dir, filename), index_col=None, header=0)
            files.append(df)

        sample_df = pd.concat(files, axis=0, ignore_index=True)

        return sample_df

    def read_all_csv_downloads(self, save=True):
        '''
        This function reads all available csv downloaded and merges into a single table. It will be used to as the pivot
        to know about all the downloaded images (and, optionally, some related metadata)

        :param save: boolean attribute whether to save the merged table into a single csv.
        :return: pd.Dataframe containing a single table for all downloads.
        '''
        print(' o Reading csv files.')
        if self.load_df and exists(self.sample_df_filepath):
            sample_df = pd.read_csv(self.sample_df_filepath)

        else:
            files = []
            for filename in listdir(self.csvs_dir):
                if filename[0] == '.': continue
                df = pd.read_csv(join(self.csvs_dir, filename), index_col=None, header=0)
                files.append(df)
            sample_df = pd.concat(files, axis=0, ignore_index=True)
            if save: sample_df.to_csv(self.sample_df_filepath)

        return sample_df

    def filter_modalities(self, sample_df, dict_file_names=None):
        '''
        This function filters the desired images as specified by containing the regular expressions found in the
        dict_file_names keys. From all the "Image Data ID" available in sample_df, it only keeps the targeted ones.

        :param sample_df: pandas dataframe containing the information from download CSV.
        :param dict_file_names: dictionary containing the BIDS entities for different filenames in the original ADNI.
        :return: pd.Dataframe -- the updated sample_df
        '''

        print(' o Reading imaging files: ', end='\r', flush=True)
        if dict_file_names is None:
            dict_file_names = ADNI_DICT_FILENAMES

        # initialise a table to be filled only with the download files we want (excluding other image modalities)
        columns = np.unique(sample_df.columns.tolist() + self.extra_download_info)
        rf = pd.DataFrame(columns=columns)

        # filter only valid rows
        valid_rows = [any([vm.lower() in adf['Description'].lower()
                           for vm in dict_file_names.keys()])
                      for _, adf in sample_df.iterrows()
                      ]

        # what has been excluded for debugging purposes.
        excluded_modalities = np.unique(sample_df[[not v for v in valid_rows]]['Description'].to_list()).tolist()

        # iterating through all images we want to keep.
        for it_image, (file_loc, image_df) in enumerate(sample_df[valid_rows].iterrows()):
            if np.mod(it_image, sum(valid_rows)//10) == 0:
                print(' o Reading imaging files: ' + str(np.round(it_image/sum(valid_rows)*100)) + '%', end='\r')

            imageid = image_df['Image Data ID']

            # BIDS ID
            subject_id = 'ADNI' + image_df['Subject'].replace('_', '')

            # Building session_id. Look at different files depending on the modality.
            # If they are not found in the corresponding file, look into ADNIMERGE
            # If still not available, it should be an error format and we provide three solutions.
            if self.pet_df is not None and (len(self.pet_df[self.pet_df['Image ID'] == int(imageid[1:])]['Visit_Month_Number']) == 1):
                adni_subject = self.adnimerge_df.loc[self.adnimerge_df['PTID'] == image_df['Subject']]
                if len(adni_subject) > 0:
                    date_format = "%Y-%m-%d"
                    scan_date = datetime.strptime(self.pet_df.loc[self.pet_df['Image ID'] == int(imageid[1:])]['Scan Date'].iloc[0], date_format)
                    diff_days = [(scan_date - datetime.strptime(d, date_format)).days for d in adni_subject['EXAMDATE']]
                    idx_adni_subject = np.argmin(np.abs(np.array(diff_days)))
                    session_id = 'm' + str(adni_subject.iloc[idx_adni_subject]['Month']).zfill(3)
                else:
                    session_id = 'm' + str(self.pet_df[self.pet_df['Image ID'] == int(imageid[1:])].iloc[0]['Visit_Month_Number']).zfill(3)

            elif self.mri_df is not None and (len(self.mri_df[self.mri_df['IMAGEUID'] == int(imageid[1:])]['Visit_Month_Number']) == 1):
                adni_subject = self.adnimerge_df.loc[self.adnimerge_df['PTID'] == image_df['Subject']]
                if len(adni_subject) > 0:
                    date_format = "%Y-%m-%d"
                    scan_date = datetime.strptime(self.mri_df.loc[self.mri_df['IMAGEUID'] == int(imageid[1:])]['SCANDATE'].iloc[0], date_format)
                    diff_days = [(scan_date - datetime.strptime(d, date_format)).days for d in adni_subject['EXAMDATE']]
                    idx_adni_subject = np.argmin(np.abs(np.array(diff_days)))
                    session_id = 'm' + str(adni_subject.iloc[idx_adni_subject]['Month']).zfill(3)
                else:
                    session_id = 'm' + str(self.mri_df[self.mri_df['IMAGEUID'] == int(imageid[1:])].iloc[0]['Visit_Month_Number']).zfill(3)

            elif self.adnimerge_df is not None and (len(self.adnimerge_df[self.adnimerge_df['IMAGEUID'] == int(imageid[1:])]['Month']) == 1):
                session_id = 'm' + str(self.adnimerge_df.loc[self.adnimerge_df['IMAGEUID'] == int(imageid[1:])].iloc[0]['Month']).zfill(3)

            else:
                sid = image_df['Visit']
                if sid in ['bl', 'sc']:
                    sid = '000'
                elif 'y' in sid:
                    year = int(sid.split('y')[1])
                    sid = str(year*12).zfill(3)
                elif isinstance(sid, int):
                    sid = str(int(sid)).zfill(3)

                session_id = 'm' + str(sid)

            # needed features for bids path construction
            file_dict = [{e: n for e,n in d.items()} for k, d in dict_file_names.items() if k.lower() in image_df['Description'].lower()][0]

            # extra info required by the program.
            image_df['subject_id'] = subject_id
            image_df['session_id'] = session_id
            image_df['dicom_path'] = imageid
            image_df['nifti_path'] = self.bids_db.build_path({**{'subject': subject_id, 'session': session_id, 'extension': 'nii.gz'}, **file_dict})
            image_df['json_path'] = self.bids_db.build_path({**{'subject': subject_id, 'session': session_id, 'extension': 'json'}, **file_dict})
            image_df['file_dict'] = file_dict

            # add to the remaining_files table
            rf = pd.concat([rf, image_df.to_frame().T])

        print(' o Reading imaging files: 100%')
        return rf, excluded_modalities

    def remove_duplicates(self, sample_df, bids_db=None):
        '''
        Rename sessions with multiple images (i.e. different Image Data ID in the same session) using the "run" entity.
        :param sample_df: pandas dataframe containing the information from download CSV (possibly filtered by the
                          desired modalities). Typically, the output of "read_all_csv_downloads" or a modified version.
        :param bids_db: (optional, use the self.bids_db as default) bids object containing the current version of the
                        dataset. Used to build the paths.
        :return: pd.Dataframe -- the updated sample_df
        '''

        if bids_db is None:
            bids_db = self.bids_db

        # DICOMS are unique, so define them as the key
        nifti_paths = {image_df['dicom_path']: image_df['nifti_path'] for _, image_df in sample_df.iterrows()}

        # Group all dicoms available for the same NIFTI filename (dict key)
        rep_im_d = {}
        for d, n in nifti_paths.items():
            if exists(n):
                # Rename existing file and add it in the first place in the list (for counting the run parameter)
                dlist = [] if n not in rep_im_d.keys() else rep_im_d[n]
                rep_im_d[n] = [d] + dlist
                t = sample_df.loc[sample_df['nifti_path'] == n].iloc[0]
                file_dict = sample_df.loc[t.name]['file_dict']
                rf = bids_db.build_path({**{
                    'subject': sample_df.loc[t.name]['subject_id'],
                    'session': sample_df.loc[t.name]['session_id'],
                    'extension': 'nii.gz',
                    'run': str(1).zfill(2)}, **file_dict})
                shutil.move(n, rf)
                shutil.move(n.replace('nii.gz', 'json'), rf.replace('nii.gz', 'json'))


            elif exists(dirname(n)) and any([all([e in n for e in f.split('_') if 'run' not in e])
                                             for f in list(filter(lambda x: 'nii.gz' in x, listdir(dirname(n))))]):

                # More than 1 file found with the nifti+run path. Adding to the list for counting in the run parameter.
                dlist = [] if n not in rep_im_d.keys() else rep_im_d[n]
                rep_im_d[n] = [f for f in list(filter(lambda x: 'nii.gz' in x, listdir(dirname(n))))
                               if all([e in n for e in f.split('_') if 'run' not in e])] + dlist
            else:
                if n in rep_im_d.keys():
                    rep_im_d[n].append(d)
                else:
                    rep_im_d[n] = [d]

        #Rename images
        rep_im_nifti = {n: dlist for n, dlist in rep_im_d.items() if len(dlist) > 1}
        for n, dlist in rep_im_nifti.items():
            total_images = sample_df.loc[sample_df['nifti_path'] == n]
            for it_t, t in enumerate(total_images.iloc):
                file_dict = sample_df.loc[t.name]['file_dict']
                sample_df.loc[t.name]['nifti_path'] = bids_db.build_path({**{
                    'subject': sample_df.loc[t.name]['subject_id'],
                    'session': sample_df.loc[t.name]['session_id'],
                    'extension': 'nii.gz',
                    'run': str(it_t+1).zfill(2)}, **file_dict})

                sample_df.loc[t.name]['json_path'] = bids_db.build_path({**{
                    'subject': sample_df.loc[t.name]['subject_id'],
                    'session': sample_df.loc[t.name]['session_id'],
                    'extension': 'json',
                    'run': str(it_t+1).zfill(2)}, **file_dict})

        return sample_df

    def remove_existing_files(self, sample_df, bids_db=None):
        '''
        Remove from sample_df the images that are already present in the dataset.
        :param sample_df: pandas dataframe containing the information from download CSV (possibly filtered by the
                          desired modalities and without duplicates).
        :param bids_db: (optional, use the self.bids_db as default) bids object containing the current version of the
                        dataset.
        :return: pd.Dataframe -- the updated sample_df
        '''

        # TODO: if multiple runs of each session (test-retest) are included in the dataset at different times,
        #  the nifti_path will not be well-specified. I.e., it will overwrite the previous image.

        if bids_db is None:
            bids_db = self.bids_db

        if 'ImageID' not in bids_db.get_entities():
            return sample_df

        # existing_image_id = [iid for iid in sample_filter_df['Image Data ID'] if len(mybids.get(ImageID=iid))>0]
        image_files = bids_db.get(extension='nii.gz')
        iid_files = [im.entities['ImageID'] for im in image_files]
        nonexisting_image_id = [iid for iid in sample_df['Image Data ID'] if iid not in iid_files]
        sample_final_df = sample_df[sample_df['Image Data ID'].isin(nonexisting_image_id)]

        # nifti_paths = {image_df['dicom_path']: image_df['nifti_path'] for _, image_df in sample_final_df.iterrows()}
        # rep_im_d = {}
        # for d, n in nifti_paths.items():
        #     if n in rep_im_d.keys():
        #         rep_im_d[n].append(d)
        #     else:
        #         rep_im_d[n] = [d]
        #
        # rep_im_nifti = {n: dlist for n, dlist in rep_im_d.items() if len(dlist) > 1}
        # for n, dlist in rep_im_nifti.items():
        #     total_images = sample_final_df.loc[sample_final_df['nifti_path'] == n]
        #     for it_t, t in enumerate(total_images.iloc):
        #         file_dict = sample_final_df.loc[t.name]['file_dict']
        #         sample_final_df.loc[t.name]['nifti_path'] = bids_db.build_path({**{
        #             'subject': sample_final_df.loc[t.name]['subject_id'],
        #             'session': sample_final_df.loc[t.name]['session_id'],
        #             'extension': 'nii.gz',
        #             'run': str(it_t + 1).zfill(2)}, **file_dict})
        #
        #         sample_final_df.loc[t.name]['json_path'] = bids_db.build_path({**{
        #             'subject': sample_final_df.loc[t.name]['subject_id'],
        #             'session': sample_final_df.loc[t.name]['session_id'],
        #             'extension': 'json',
        #             'run': str(it_t + 1).zfill(2)}, **file_dict})

        return sample_final_df

    def unzip_dicom(self, sample_df, zipfiles=None, bids_db=None):
        '''

        :param sample_df: pandas dataframe containing the information from download CSV (possibly filtered by the
                          desired modalities and without duplicates).
        :param zipfiles: list of zips to be decompressed
        :param bids_db: (optional, use the self.bids_db as default) bids object containing the current version of the
                        dataset.
        :return: pd.Dataframe, an updated version of sample_df.
        '''

        unique_iids = sample_df['Image Data ID'].to_numpy()
        new_sample = pd.DataFrame(columns=sample_df.columns.tolist())


        if zipfiles is None:
            zipfiles = self.zips_download

        if not isinstance(zipfiles, list):
            zipfiles = [zipfiles]

        zipfile_dict = {}

        # Iterate over all zipfiles
        for filename in zipfiles:
            if filename[0] == '.': continue

            # Read zip file
            zip = ZipFile(join(self.zips_dir, filename))

            # Get filepaths to each image folder (i.e., contiaining the zips). Zip.namelist provide the path to each
            # DICOM image, but we want only its directory
            unique_filepaths_zip = np.unique([dirname(z) for z in zip.namelist()])

            # Filter those images to decompress within the zip file (i.e., that are in the sample_df dataframe)
            files = {basename(fp): fp for fp in unique_filepaths_zip if basename(fp) in unique_iids}
            zipfile_dict[filename] = {'zipname': filename, 'zip': zip, 'files': files}
            for it_file, (iid, fp) in enumerate(files.items()):
                if np.mod(int(it_file/len(files)*100), 10) == 0:
                    print('  o Unzipped files: ' + str(int(it_file/len(files)*100)) + '%', end='\r')

                # Get IMAGE Instance
                image_df = sample_df.loc[sample_df['Image Data ID'] == iid]# if iid not in new_sample.index: continue
                if len(image_df) > 0:
                    image_df = image_df.iloc[0]
                else:
                    continue

                # Unzip if it doesn't exist in the WC directory
                image_df['dicom_path'] = join(self.working_dir, fp)
                if not exists(join(self.working_dir, fp)):
                    makedirs(join(self.working_dir, fp))
                    subprocess.run(['unzip', join(self.zips_dir, filename), join(fp, '*'), '-d', join(self.working_dir, fp)],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                # Update the nifti_path and json_path for PET if the tracer is not specified
                if 'pet' in image_df['nifti_path'] and 'trc' not in image_df['nifti_path']:
                    subject_id = 'ADNI' + replace_underscore(image_df['Subject'])
                    session_id = basename(image_df['nifti_path']).split('_')[1].split('-')[1]
                    file_dict = [{e: n for e, n in d.items()}
                                 for k, d in ADNI_DICT_FILENAMES.items()
                                 if k.lower() in image_df['Description'].lower()][0]

                    image_df['nifti_path'] = bids_db.build_path({**{
                        'subject': subject_id,
                        'session': session_id,
                        'extension': 'nii.gz'}, **file_dict})

                    image_df['json_path'] = bids_db.build_path({**{
                        'subject': subject_id,
                        'session': session_id,
                        'extension': 'json'}, **file_dict})

                new_sample = pd.concat([new_sample, image_df.to_frame().T])

        print('  o Unzipped files: 100%')
        return new_sample


class DCM2NII(object):

    def __init__(self, mri_df, pet_df):
        self.mri_df = mri_df
        self.pet_df = pet_df

    def _dcm2niix(self, image_df, num_processed=None, num_total=None):
        '''
        This function calls the dcm2niix function. First, it uses the FSL by default. If it does not work, it uses
        the github function from rodenlab.

        :param image_df: (pd.Dataframe) contains the information from download CSV (possibly filtered by the
                          desired modalities and without duplicates).
        :param num_processed: (int, optional). Number of already processed images, for tracking (i.e. print) purposes
                              Default=0
        :param num_total: (int, optional). Number of images to be converted, for tracking (i.e. print) purposes.
                          Default=len(sample_df):param num_processed:
        :return: failed_images (dict) of images that could not be converted with the error message.
        '''

        if np.mod(num_processed, 5) == 0:
            print('* Converting images DCM to NIFTI: ' + str(num_processed) + '/' + str(num_total), end='\r')

        dicom_path = image_df['dicom_path']
        nifti_path = image_df['nifti_path']
        json_path = image_df['json_path']

        if not exists(dirname(nifti_path)):
            makedirs(dirname(nifti_path))

        if isdir(dicom_path):
            dcmout = subprocess.run(["dcm2niix", "-z", "y", "-o", dirname(nifti_path), dicom_path],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            dcmerr = dcmout.stderr.decode("utf-8")
            dcmcode = dcmout.returncode
            if dcmcode != 0 or dcmerr != '' or 'No valid' in dcmout.stdout.decode("utf-8"):
                if exists(DCM2NIIX_GITHUB):
                    dcmout = subprocess.run(
                        [DCM2NIIX_GITHUB, "-b", "y", "-w", "0", "-z", "y", "-o", dirname(nifti_path), dicom_path],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    dcmerr = dcmout.stderr.decode("utf-8")
                    dcmcode = dcmout.returncode
                    if dcmerr != '':
                        failed_images = {'dicom_path': dicom_path, 'error': dcmerr, 'code': dcmcode}
                        return failed_images

                    elif 'No valid' in dcmout.stdout.decode("utf-8"):
                        failed_images = {'dicom_path': dicom_path, 'error': "No valid DICOM files", 'code': dcmcode}
                        return failed_images

                    elif dcmcode != 0:
                        failed_images = {'dicom_path': dicom_path, 'error': dcmout.stdout.decode("utf-8"), 'code': dcmcode}
                        return failed_images


                else:
                    failed_images = {'dicom_path': dicom_path, 'error': dcmerr, 'code': dcmcode}
                    return failed_images

            # Rename .nii.gz and .json files
            self.rename_image(nifti_path, json_path, image_df['Image Data ID'])

            # Add relevant info to .json file
            self.fill_json(image_df, json_path)

            return None

        else:

            failed_images = {'dicom_path': dicom_path, 'error': 'zip has not decompressed well.'
                                                                'Try to decompress it manually', 'code': 'retry'}

            return failed_images

    def convert(self, sample_df, num_processed=None, num_total=None, num_cores=1):
        '''
        This function tries to convert DICOM images to NIFTI standard. Sometimes it fails, for a variety of reasons
        specified in the fields of dictionary "failed_images".

        :param sample_df: (pd.Dataframe) contains the information from download CSV (possibly filtered by the
                          desired modalities and without duplicates).
        :param num_processed: (int, optional). Number of already processed images, for tracking (i.e. print) purposes
                              Default=0
        :param num_total: (int, optional). Number of images to be converted, for tracking (i.e. print) purposes.
                          Default=len(sample_df)
        :param num_cores: (int, optional). Number of cores for parallel processing. Default=1 (no parallel processing)
        :return: failed_images (list of dicts) images that could not be converted with the error message.
        '''

        if num_total is None:
            num_total = len(sample_df)
        if num_processed is None:
            num_processed = 0

        if num_cores > 1:
            results = Parallel(n_jobs=num_cores)(
                delayed(self._dcm2niix)(image_df, num_processed=it_df, num_total=num_total)
                for it_df, (_, image_df) in enumerate(sample_df.iterrows())
            )
            failed_images = [r for r in results if r is not None]

        else:
            failed_images = []
            for it_df, (iid, image_df) in enumerate(sample_df.iterrows()):
                if 'Image Data ID' not in image_df:
                    image_df['Image Data ID'] = iid

                fim = self._dcm2niix(image_df, num_processed=num_processed+it_df, num_total=num_total)
                if fim is not None: failed_images += [fim]

        return failed_images

    def rename_image(self, nifti_filepath, json_filepath, image_id):
        '''
        This function renames images --nifti and json-- (i.e. moves if different directory) from the converted directory to the
        new directory and named accordingly.

        :param nifti_filepath: str filepath specifying the .nii.gz
        :param json_filepath: str filepath specifying the .json
        :param image_id: image identifier as specified in ADNI
        :return: None
        '''

        images_in_folder = listdir(dirname(nifti_filepath))
        for image in images_in_folder:
            if image_id not in str(image):
                continue

            elif 'json' in image:
                rename(join(dirname(nifti_filepath), image), json_filepath)

            elif 'nii.gz' in image:
                rename(join(dirname(nifti_filepath), image), nifti_filepath)

            else:
                print('[Warning] There is a non .json and non .nii.gz file in the following directory: ' + nifti_filepath)

    def fill_json(self, image_df, json_filepath):
        '''
        This function add new information to the JSON file, such as the Imge ID.

        :param image_df: pandas dataframe containing metadata from the image associated to this json filepath.
        :param json_filepath: str specifying the json filepath.
        :return: None
        '''
        json_file = open(json_filepath, 'r')
        json_dict = json.load(json_file)
        json_dict['ImageID'] = image_df['Image Data ID']

        if image_df['Modality'].lower() == 'pet':
            json_dict['TracerName'] = [d['tracer'] for k, d in ADNI_DICT_FILENAMES.items() if
                                       k.lower() in image_df['Description'].lower() and 'tracer' in d][0]
            json_dict['TracerRadionuclide'] = 'TracerRadionuclide'
            json_dict['ModeOfAdministration'] = 'ModeOfAdministration'
            json_dict['InjectedRadioactivity'] = 666
            json_dict['InjectedRadioactivityUnits'] = 'MBq'
            json_dict['InjectedMass'] = 666
            json_dict['InjectedMassUnits'] = 'ug'
            json_dict['SpecificRadioactivity'] = 666
            json_dict['SpecificRadioactivityUnits'] = 'MBq/ug'
            json_dict['TimeZero'] = '00:00:00'
            json_dict['ScanStart'] = 0
            json_dict['InjectionStart'] = 0
            json_dict['FrameTimesStart'] = '' if 'FrameTimesStart' not in json_dict else json_dict['FrameTimesStart']
            json_dict['AcquisitionMode'] = 'list mode'
            json_dict['ImageDecayCorrected'] = True
            json_dict['ImageDecayCorrectionTime'] = 0
            json_dict['ReconMethodName'] = '"3D-OSEM-PSF'
            json_dict['ReconMethodParameterLabels'] = ["subsets", "iterations"]
            json_dict['ReconMethodParameterUnits'] = ["none", "none"]
            json_dict['ReconMethodParameterValues'] = [666, 666]
            json_dict['ReconFilterType'] = 'none'
            json_dict['ReconFilterSize'] = 666
            json_dict['AttenuationCorrection'] = "[137Cs]transmission scan-based"
            json_dict['InfusionRadioactivity'] = 666
            json_dict['InfusionStart'] = 666
            json_dict['InfusionSpeed'] = 666
            json_dict['InfusionSpeedUnits'] = 'InfusionSpeedUnits'
            json_dict['InjectedVolume'] = 666
            json_dict['ScanDate'] = ''
            if len(self.pet_df.loc[self.pet_df['Image ID'] == int(json_dict['ImageID'][1:])]) > 0:
                idx = self.pet_df['Image ID'] == int(json_dict['ImageID'][1:])
                json_dict['ScanDate'] = self.pet_df.loc[idx].iloc[0]['Scan Date']

            json_object = json.dumps(json_dict, indent=4)
            json_file = open(json_filepath, "w")
            json_file.write(json_object)

        elif image_df['Modality'].lower() == 'fmri':
            json_dict['TaskName'] = 'Resting State'
            json_dict['ScanDate'] = ''
            if len(self.mri_df.loc[self.mri_df['IMAGEUID'] == int(json_dict['ImageID'][1:])]) > 0:
                idx = self.mri_df['IMAGEUID'] == int(json_dict['ImageID'][1:])
                json_dict['ScanDate'] = self.mri_df.loc[idx].iloc[0]['SCANDATE']

            json_object = json.dumps(json_dict, indent=4)
            json_file = open(json_filepath, "w")
            json_file.write(json_object)

        elif image_df['Modality'].lower() == 'mri':
            json_dict['ScanDate'] = ''
            if len(self.mri_df.loc[self.mri_df['IMAGEUID'] == int(json_dict['ImageID'][1:])]) > 0:
                idx = self.mri_df['IMAGEUID'] == int(json_dict['ImageID'][1:])
                json_dict['ScanDate'] = self.mri_df.loc[idx].iloc[0]['SCANDATE']

            json_object = json.dumps(json_dict, indent=4)
            json_file = open(json_filepath, "w")
            json_file.write(json_object)

        else:
            print('Image ' + image_df['Image Data ID'] + ' unwanted image modality (' +
                  image_df['Modality'] + ') in folder. Subject:' + image_df['subject_id'])


def create_dataset_description(main_folder_path, dataset_name='ADNI'):
    '''
    This function creates a .json file with the following fields:
    Name, BIDSVersion, License, Authors, Acknowledgements, HowToAcknowledge,
    Funding, ReferencesAndLinks and DatasetDOI

    :param main_folder_path:  where to store the dataset_description.json file.
    :param dataset_name: str naming the dataset.
    :return: None
    '''

    data_descr = {}
    data_descr['Name'] = dataset_name
    data_descr['BIDSVersion'] = ''
    data_descr['License'] = ''
    data_descr['Authors'] = ['']
    data_descr['Acknowledgements'] = ''
    data_descr['HowToAcknowledge'] = ''
    data_descr['Funding'] = ['']
    data_descr['ReferencesAndLinks'] = ['']
    data_descr['DatasetDOI'] = ''
    data_descr_path = join(main_folder_path, 'dataset_description.json')

    json_object = json.dumps(data_descr, indent=4)
    with open(data_descr_path, 'w') as outfile:
        outfile.write(json_object)

def init_participants(sample_df, bids_db, adnimerge_df, specific_cols=None):
    '''
    This functions just calls two different functions to create the participants tsv and json files.
    Before quitting, it also creates subjects folders.

    :param sample_df: pandas dataframe containing the information from download CSV (possibly filtered by the
                      desired modalities and without duplicates).
    :param bids_db: (optional, use the self.bids_db as default) bids object containing the current version of the
                    dataset.
    :param adnimerge_df: pandas dataframe containing adnimerge file
    :param specific_cols: (optional) what columns from ADNIMERGE to save. If None, we will add
                          ['sex', 'ethnic', 'race', 'apoe4', 'education', 'group_bl', 'age_bl', 'Married']
    :return:
    '''
    '''

    This function calls the two functions 'create_participants_tsv' and
    'create_participants_json' that create the .tsv and .json participants files

    Parameters
    ----------
    main_folder_path : TYPE str
        DESCRIPTION. Path to the main folder where all BIDS structure is contained.
    info_data : TYPE str
        DESCRIPTION. Path and name of the .csv file containing the downloaded ADNI dataset metadata.

    Returns
    -------
    subjs : TYPE NumPy Object array
        DESCRIPTION. List containing the participant Id of each subject in the .csv file pointed in
        the variable 'info_data' Example: 003_S_6264
        This variable is created in 'crear_participants_tsv' function and is used in:
            subjects Folder, subject var (subject_id_sessions.tsv)
            and folder cration (Find sessions that we have images from to only create those)
    info_csv : TYPE DataFrame
        DESCRIPTION. Matrix containing all metadata of the downloaded ADNI dataset.
        It's columns are 'Data ID'',"Subject","Group","Sex","Age","Visit","Modality",
        "Description","Type","Acq Date","Format","Downloaded"
        This variable is created in 'crear_participants_tsv' function and is used in:
            Folder cration (Find sessions that we have images from to only create those)

    '''
    participants_cols = ['participant_id', 'num_subject', 'adni_id']

    if specific_cols is None:
        specific_cols = PARTICIPANTS_COL

    init_participants_json(bids_db.root, participants_cols + specific_cols)
    pdict = init_participants_tsv(sample_df, bids_db, adnimerge_df, participants_cols, specific_cols)

    for sid in pdict.keys():
        if not isdir(join(bids_db.root, 'sub-' + sid)):
            makedirs(join(bids_db.root, 'sub-' + sid))

    return pdict

def init_participants_tsv(sample_df, bids_db, adnimerge_df, participants_cols, specific_cols):
    '''
    This function initializes the participants tsv file. It first looks if exist: if it does, it reads everything and
    keep adding new unexisting information (in practice it overwrites); if it doesn't, it directly writes everything
    new.

    :param sample_df: pandas dataframe containing the information from download CSV (possibly filtered by the
                      desired modalities and without duplicates).
    :param bids_db: (optional, use the self.bids_db as default) bids object containing the current version of the
                    dataset.
    :param adnimerge_df: pandas dataframe containing adnimerge file
    :param participants_cols: required participants cols
    :param specific_cols: optional participants cols, more related to each dataset
    :return: tuple(dict, dict): the first dictionary is the existing subjects in the dataset, the second is the missing
                                subjects in the dataset (to be written extra).
    '''
    '''

    This function creates a .tsv file with the following column headers:
        'participant_id','SubjectID', 'AlternativeID', 'Group', 'Sex' and 'Age'

    Parameters
    ----------
    info_data : TYPE str
        DESCRIPTION. Path and name of the .csv file containing the downloaded ADNI dataset metadata.
    main_folder_path : TYPE str
        DESCRIPTION. Path to the main folder where all BIDS structure is contained.

    Returns
    -------
    subjs : TYPE NumPy Object array
        DESCRIPTION. List containing the participant Id of each subject in the .csv file pointed in
        the variable 'info_data' Example: 003_S_6264
        This variable is created in 'crear_participants_tsv' function and is used in:
            subjects Folder, subject var (subject_id_sessions.tsv)
            and folder cration (Find sessions that we have images from to only create those)
    info_csv : TYPE DataFrame
        DESCRIPTION. Matrix containing all metadata of the downloaded ADNI dataset.
        It's columns are 'Data ID'',"Subject","Group","Sex","Age","Visit","Modality",
        "Description","Type","Acq Date","Format","Downloaded"
        This variable is created in 'crear_participants_tsv' function and is used in:
            Folder cration (Find sessions that we have images from to only create those)

    '''

    subject_list = sample_df['Subject'].unique()

    col_names = participants_cols + specific_cols

    existing_p_dict = {}
    if exists(join(bids_db.root, 'participants.tsv')):
        p_df = pd.read_csv(join(bids_db.root, 'participants.tsv'), delimiter='\t')
        existing_p_dict = {p['participant_id']: p.to_dict() for _, p in p_df.iterrows()}
        existing_subjects = [p['adni_id'] for pid, p in existing_p_dict.items()]
        subject_list = existing_subjects + [s for s in subject_list if s not in existing_subjects]


    to_write = list(existing_p_dict.values())
    participants_dict = {}
    for it_sid, sid in enumerate(subject_list):
        sid_bids = 'ADNI' + replace_underscore(sid)

        if sid_bids not in existing_p_dict.keys():
            d = {
                'participant_id': sid_bids,
                'num_subject': str(it_sid).rjust(5, '0'),
                'adni_id': sid,
            }

            adni_subject = adnimerge_df.loc[adnimerge_df['PTID'] == sid]
            if len(adni_subject) > 0:
                baseline_visit = np.argmin(adnimerge_df.loc[adnimerge_df['PTID'] == sid]['AGE'])
                metadata_bl = adnimerge_df.iloc[baseline_visit]
                d = {**d, **{bids_k: metadata_bl[ADNI_MERGE_DICT[bids_k]] for bids_k in specific_cols}}

            participants_dict[sid_bids] = d


            to_write.append(d)

    if to_write:
        with open(join(bids_db.root, 'participants.tsv'), 'w') as csv_file:
            tsv_writer = csv.DictWriter(csv_file, fieldnames=col_names, delimiter='\t')
            tsv_writer.writeheader()
            tsv_writer.writerows(to_write)

    return {**existing_p_dict, **participants_dict}

def init_participants_json(bids_dir, col_names):
    '''
    This function creates a .json file explaining the data fields you can find to the participants tsv file.

    :param bids_db: (optional, use the self.bids_db as default) bids object containing the current version of the
                    dataset.
    :param col_names:
    :return: None
    '''

    if all([c in list(participants_description.keys()) for c in col_names]):
        data_descr_path = join(bids_dir, 'participants.json')
        json_object = json.dumps(participants_description, indent=4)
        with open(data_descr_path, 'w') as outfile:
            outfile.write(json_object)
    else:
        raise ValueError('Participants JSON file missing some columns. Please, go to the code and re-write the '
                         'specifications before running the code again. Missing columns: '
                         + str([c for c in col_names if c not in list(participants_description.keys())]))

def init_sessions(p_dict, sample_df, bids_db, adnimerge_df, specific_cols=None):
    '''

    :param p_dict: pd.Dataframe with the participants table.
    :param sample_df: pandas dataframe containing the information from download CSV (possibly filtered by the
                      desired modalities and without duplicates).
    :param bids_db: (optional, use the self.bids_db as default) bids object containing the current version of the
                    dataset.
    :param adnimerge_df: pandas dataframe containing adnimerge file
    :param specific_cols: columns from the ADNIMERGE file that we want to store at the session level.
    :return: None
    '''
    '''
    This function creates .tsv file that contains the metadata of all sessions

    Parameters
    ----------
    subjs : TYPE NumPy Object array
        DESCRIPTION. List containing the participant Id of each subject in the .csv file pointed in
        the variable 'info_data' Example: 003_S_6264
        This variable is created in 'crear_participants_tsv' function and is used in:
            subjects Folder, subject var (subject_id_sessions.tsv)
            and folder cration (Find sessions that we have images from to only create those)
    main_folder_path : TYPE str
       DESCRIPTION. Path to the main folder where all BIDS structure is contained.
    subjects_prefix : TYPE str
        DESCRIPTION.Text prefix with which each subject file has been named 'sub-'
    visits_df : TYPE dataframe
        DESCRIPTION. Dictionary containing the visits in the ADNI protocol and their code.
        Columns: 	Phase	ID	VISCODE	VISNAME	VISORDER
    idas_df : TYPE dataframe
        DESCRIPTION. Dataframe containing the information of the 'idasearch.csv' file, downloaded from ida webpage
        It contains the information on all subjects, subject's sessions and their description
        Columns: Subject_ID, Projet, Phase, Sex, Weight, Research_Group, APOE A1, APOE A2, Visit, Study Date, Archive Date, Age, Modality, Description, Imaging Protocol.
    mrilist_df : TYPE dataframe
        DESCRIPTION. List of all MRI images and their description.
        Columns: Unnamed: 0(Index), TYPE, SUBJECT, VISIT, MAGSTRENGTH, SEQUENCE, SCANDATE, STUDYID, IMAGEUID, Visit_Month_Number
    petlist_df : TYPE dataframe
       DESCRIPTION. List of all PET images and their description.
       Columns: Unnamed: 0(Index), TYPE, SUBJECT, VISIT, MAGSTRENGTH, SEQUENCE, SCANDATE, STUDYID, IMAGEUID, Visit_Month_Number

    Returns
    -------
    None.

    '''

    session_cols = ['session_id', 'session_description']
    if specific_cols is None: specific_cols = SESSION_COLS


    col_names = session_cols + specific_cols
    sess_dict = {}

    # Iterate over all participants
    for it_s, (sid, s_dict) in enumerate(p_dict.items()):

        sess_desc_path_json = join(bids_db.root, 'sub-' + sid, 'sub-' + sid + '_sessions.json')
        sess_desc_path_tsv = join(bids_db.root, 'sub-' + sid, 'sub-' + sid + '_sessions.tsv')

        json_object = json.dumps(session_description, indent=4)
        with open(sess_desc_path_json, 'w') as outfile:
            outfile.write(json_object)

        sess_dict[sid] = {}
        # if the subject is not in the final dataframe with the images to convert, skip it.
        if sid not in sample_df['subject_id'].unique():
            to_write = []

        else:
            subject_df = sample_df[sample_df['Subject'] == s_dict['adni_id']]
            session_list = subject_df['session_id'].unique()

            # if already exists (in such case means we are adding a new session to an existing participant), read it.
            sessions_df = bids_db.get(subject=sid, extension='tsv', regex_search='sessions')
            if sessions_df:
                sess_dict[sid] = {sinfo['session_id']: sinfo for sinfo in bids_db.get_file(sessions_df[0].path).get_df().to_dict(orient='records')}

            to_write = list(sess_dict[sid].values())

            # iterate over all sessions available for the participant and added to the file
            for sess_viscode in session_list:
                sess_description = subject_df['Visit'].unique()
                session_adnimerge = adnimerge_df.loc[(adnimerge_df['PTID'] == s_dict['adni_id']) &
                                                     (adnimerge_df['Month'] == int(sess_viscode.split('m')[1]))]

                d = {
                    'session_id': sess_viscode,
                    'session_description': '_'.join(sess_description)
                }

                if len(session_adnimerge) > 0:
                    sess_metadata = session_adnimerge.iloc[0]

                    d_s = {s: sess_metadata[ADNI_MERGE_DICT[s]] for s in specific_cols if s in ADNI_MERGE_DICT.keys()}
                    d_s['time_to_bl_days'] = sess_metadata['Years_bl'] * 365.25
                    d_s['age'] = sess_metadata['AGE'] + sess_metadata['Years_bl']

                else:
                    d_s = {s: '' for s in specific_cols}

                to_write.append({**d, **d_s})

        if to_write:
            with open(sess_desc_path_tsv, 'w') as csv_file:
                tsv_writer = csv.DictWriter(csv_file, fieldnames=col_names, delimiter='\t')
                tsv_writer.writeheader()
                tsv_writer.writerows(to_write)
