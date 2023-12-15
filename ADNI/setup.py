from os.path import join, dirname

###################
# Input variables #
###################
path_data = '/media/acasamitjana/BIG_DATA/'
path_utils = '/home/acasamitjana/data/util-files'
bids_dir = '/media/acasamitjana/BIG_DATA/ProvaRepo/ADNI-BIDS/rawdata'
db_file = join(dirname(bids_dir), 'Prova-Repo.db')
db_file_wc = join(dirname(bids_dir), 'Prova-Repo-wc.db')

DCM2NIIX_GITHUB = '/home/acasamitjana/software/dcm2niix/build/bin/dcm2niix'


##########################
# ADNI related variables #
##########################
ADNI_DICT_FILENAMES = {
    # 'PET': 'pet',
    'fMRI': {'suffix': 'bold', 'task': 'rest'},
    'fcMRI': {'suffix': 'bold', 'task': 'rest'},
    # 'Hipp': 't2hipp',
    'GRAPPA': {'acquisition': 'acc', 'suffix': 'T1w'},
    'SENSE': {'acquisition': 'acc', 'suffix': 'T1w'},
    'Accelerated_Sag': {'acquisition': 'acc', 'suffix': 'T1w'},
    'MPRAGE': {'suffix': 'T1w'},
    'MP-RAGE': {'suffix': 'T1w'},
    'IR-SPGR': {'suffix': 'T1w'},
    'IR-FSPGR': {'suffix': 'T1w'},
    'AV45': {'tracer': 'fbp', 'suffix': 'pet'},
    'AV-45': {'tracer': 'fbp', 'suffix': 'pet'},
    'AV_45': {'tracer': 'fbp', 'suffix': 'pet'},
    'florbetapir': {'tracer': 'fbp', 'suffix': 'pet'},
    'FBB': {'tracer': 'fbb', 'suffix': 'pet'},
    'Florbetaben': {'tracer': 'fbb', 'suffix': 'pet'},
    'PIB': {'tracer': 'pib', 'suffix': 'pet'},
    'FLAIR': {'suffix': 'FLAIR'},
    'T2': {'suffix': 'T2w'},
}

ADNI_MERGE_DICT = {
    'sex': 'PTGENDER',
    'education': 'PTEDUCAT',
    'apoe4': 'APOE4',
    'ethnic': 'PTETHCAT',
    'race': 'PTRACCAT',
    'age_bl': 'AGE',
    'age': 'AGE',
    'group_bl': 'DX',
    'time_to_bl_months': 'Month_bl',
    'time_to_bl_years': 'Years_bl',
    'dx': 'DX',
    'phase': 'COLPROT',
    'site': 'SITE',
    'viscode': 'VISCODE',
    'icv': 'ICV',
    'mmse': 'MMSE',
    'fdg': 'FDG',
    'abeta': 'ABETA',
    'ptau': 'PTAU',
    'ttau': 'TAU',
    'cdrsb': 'CDRSB',
    'adas11': 'ADAS11',
    'adas13': 'ADAS13',
    'faq': 'FAQ',
    'Married': 'PTMARRY',
    'HPVol': 'Hippocampus'
}

DX_LONG_DICT = {
    -1: 'Unknown',
    0: 'CN_cross',
    1: 'MCI_cross',
    2: 'Dementia_cross',
    3: 'CN_stable',
    4: 'CN_MCI_converter',
    5: 'CN_AD_converter',
    6: 'MCI_stable',
    7: 'MCI_AD_converter',
    8: 'AD_stable',
    9: 'Reversed_diagnosis',
}

def get_dx(dx):
    if dx == '':
        dx = -1
    elif dx == '-1':
        dx = -1
    elif dx in ['CN', 'Control']:
        dx = 0
    elif dx == 'MCI':
        dx = 1
    elif dx in ['Dementia', 'AD']:
        dx = 2

    return dx

def get_dx_long(dx_list):
    if len(dx_list) == 1:
        return get_dx(dx_list[0])

    dx_bl = get_dx(dx_list[0])
    dx_last = get_dx(dx_list[-1])

    if dx_bl == -1 or dx_last == -1:
        dx_flag = [get_dx(d) != -1 for d in dx_list]

        if sum(dx_flag) == 0:
            return -1

        elif sum(dx_flag) == 1:
            idx = int(np.where(np.array(dx_flag) == 1)[0])
            return get_dx(dx_list[idx])

        else:
            idx = np.where(np.array(dx_flag) == 1)[0]
            idx_min = int(np.min(idx))
            idx_max = int(np.min(idx))
            dx_bl = get_dx(dx_list[idx_min])
            dx_last = get_dx(dx_list[idx_max])

    if dx_bl == 0:
        if dx_last == 0:
            dx = 3  # CN stable
        elif dx_last == 1:
            dx = 4  # CN-MCI converter
        elif dx_last == 2:
            dx = 5  # CN-AD converter
        else:
            dx = 9
    elif dx_bl == 1:
        if dx_last == 1:
            dx = 6  # MCI_stable
        elif dx_last == 2:
            dx = 7  # MCI-AD converter
        else:
            dx = 9
    elif dx_bl == 2:
        if dx_last == 2:
            dx = 8  # AD stable
        else:
            dx = 9
    else:
        dx = 9  # Reversed diagnosis

    return dx

################################
# Variables used to build BIDS #
################################

PARTICIPANTS_COL = ['sex', 'ethnic', 'race', 'apoe4', 'education', 'group_bl', 'age_bl', 'Married']

participants_description = {
    "participant_id": {
        "Description": "Participant identifier in BIDS format.",
    },

    "adni_id": {
        "Description": "Identification as downloaded from ADNI Dataset"
    },

    "num_subject": {
        "Description": "Number of subject ordered by the inclusion in the BIDS dataset.",
    },

    "sex": {
        "Description": "sex of the participant as reported by the participant",
        "Levels": {
            "M": "male",
            "F": "female"
        }
    },

    "group_bl": {
        "Description": 'Baseline diagnostic label.',
        "Levels": {
            "CN": '?'  ##### WHat are the levels?
        }
    },

    "age_bl": {
        "Description": "Age of the participant at the baseline visit",
        "Units": "years"
    },

    "ethnic": {
        "Description": ""
    },

    "race": {
        "Description": ""
    },

    "apoe4": {
        "Description": "Number of apoe4 alelles",
        "Levels": {
            "0": "It can be a2a2, a2a3, a3a2, a3a3",
            "1": "It can be a2a4, a4a2, a3a4, a4a3",
            "2": "It's a4a4"
        }
    },

    "education": {
        "Description": "Years of education"
    },

    "Married": {
        "Description": "Married, never married, divorced, etc... at baseline."
    },

}

SESSION_COLS = ['age', 'time_to_bl_days', 'time_to_bl_months', 'time_to_bl_years', 'dx', 'phase', 'site',
                'viscode', 'icv', 'fdg', 'abeta', 'ptau', 'ttau', 'mmse', 'cdrsb', 'adas11', 'adas13', 'faq', 'HPVol']

session_description = {}
session_description['session_id'] = 'session identifier'
session_description['session_description'] = 'session description as in VISITS.csv file from the ADNI database'
session_description['age'] = ''
session_description['time_to_bl_days'] = ''
session_description['time_to_bl_months'] = ''
session_description['time_to_bl_years'] = ''
session_description['dx'] = ['']
session_description['phase'] = ''
session_description['site'] = ''
session_description['viscode'] = ''
session_description['icv'] = ''
session_description['fdg'] = ''
session_description['abeta'] = ''
session_description['ptau'] = ''
session_description['ttau'] = ''
session_description['mmse'] = ''
session_description['cdrsb'] = ''
session_description['adas11'] = ''
session_description['adas13'] = ''
session_description['faq'] = ''
session_description['HPVol'] = ''
