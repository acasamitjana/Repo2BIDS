import pdb
from os.path import join, dirname, exists
from os import makedirs, remove
import bids
from argparse import ArgumentParser
import shutil
import subprocess

from utils import DCM2NII, ADNI_DB, create_dataset_description, init_participants, init_sessions, final_report
from ADNI.setup import *

################
# Let's start! #
################
arg_parser = ArgumentParser(description='Computes the prediction of certain models')
arg_parser.add_argument('--c', default=None, nargs='+', help="CSV to read")
arg_parser.add_argument('--z', default=None, nargs='+', help="Zipfiles to decompress")
arg_parser.add_argument('--no_ovrw', action='store_true', help="If don't want to overwrite the existing BIDS DB")
args = arg_parser.parse_args()

# Reading BIDS dataset
if not exists(bids_dir): makedirs(bids_dir)
create_dataset_description(bids_dir)

bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
if not exists(db_file):
    print('Initialising the BIDS dataset.')
    if not exists(bids_dir): makedirs(bids_dir)
    create_dataset_description(bids_dir)

    bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
    # bids_db.save(db_file)

else:
    print('Initialising the BIDS dataset. Reading from ' + db_file)
    if exists(db_file_wc): shutil.rmtree(db_file_wc)
    shutil.move(db_file, db_file_wc)
    bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False, database_path=db_file_wc)

adni_db = ADNI_DB(path_data, path_utils, bids_db, load_df=True)
adni_db.build()

init_df = adni_db.read_all_csv_downloads(save=True) if args.c is None else adni_db.read_csvs(args.c)
sample_df, excluded_mod = adni_db.filter_modalities(sample_df=init_df)
sample_df = adni_db.remove_existing_files(sample_df=sample_df)
sample_df = adni_db.remove_duplicates(sample_df=sample_df)

failed_images = []
pdb.set_trace()
zipfiles = adni_db.zips_download if args.z is None else args.z
for it_z, zipfile in enumerate(zipfiles):

    print('\nUnzipping image files.')
    sample_it_df = adni_db.unzip_dicom(sample_df=sample_df, zipfiles=zipfile)

    print('\nWriting participants file.')
    p_dict = init_participants(sample_it_df, bids_db, adni_db.adnimerge_df)

    print('\nCreating sessions tsv files.')
    init_sessions(p_dict, sample_it_df, bids_db, adni_db.adnimerge_df)

    print('\nConverting DICOM to NIFTI images.')
    DCM = DCM2NII(adni_db.mri_df, adni_db.pet_df)
    f = DCM.convert(sample_it_df)
    failed_images.extend(f)

print('\nWriting final report')
final_report(init_df, sample_df, failed_images)


if not args.no_ovrw:
    if exists(db_file): shutil.rmtree(db_file)
    bids_db = bids.layout.BIDSLayout(root=bids_dir, validate=False)
    bids_db.save(db_file)
