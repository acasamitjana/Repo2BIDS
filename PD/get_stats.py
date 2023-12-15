import csv
import pdb
from os.path import join

import bids
import numpy as np

def get_int(a):
    try:
        return int(a)
    except:
        return a

def read_sess_file(file):
    def get_value(a):
        try:
            return float(a)
        except:
            return a

    with open(file, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter='\t')
        data = {k: [] for k in csvreader.fieldnames}
        for row in csvreader:
            data = {k: data[k] + [get_value(row[k])] for k in csvreader.fieldnames}

    return data


#https://peerherholz.github.io/workshop_weizmann/advanced/pybids.html
bids_dir = '/vol/home/student_pd12/ADRIA/PD-BIDS'

print('Reading dataset.')
mybids = bids.layout.BIDSLayout(root=bids_dir, validate=False)

stats_dict = {}
participants_file = mybids.get(return_type='filename', extension='.tsv', regex_search='participants')
participants_csv = read_sess_file(participants_file[0])

subject_list = mybids.get_subjects()
for subject in subject_list:
    # pdb.set_trace()
    idx_subject = np.where(np.array(participants_csv['participant_id']) == subject)[0][0]
    sess_file = mybids.get(return_type='filename', subject=subject, extension='.tsv', regex_search='sessions')
    sess_csv = read_sess_file(sess_file[0])
    session_list = mybids.get(subject=subject, return_type='id', target='session')
    idx_bl = np.argmin(np.array(sess_csv['age']))

    stats_dict[subject] = {}
    stats_dict[subject]['num_sessions'] = len(session_list)
    stats_dict[subject]['sex'] = participants_csv['sex'][idx_subject]
    stats_dict[subject]['education'] = get_int(participants_csv['education'][idx_subject])
    stats_dict[subject]['age_bl'] = float(min([float(a) for a in sess_csv['age']]))
    stats_dict[subject]['dx'] = participants_csv['diagnosis'][idx_subject]
    stats_dict[subject]['ages'] = [float(a) for a in sess_csv['age']]
    stats_dict[subject]['timespan'] = float(max([float(t) for t in sess_csv['time_to_bl_days']]) / 365.25)
    stats_dict[subject]['total_images'] = len(mybids.get(subject=subject, extension='nii.gz'))
    stats_dict[subject]['imaging'] = {
        't1w': sum([1 for sess in session_list if len(mybids.get(subject=subject, extension='nii.gz', suffix='T1w', session=sess)) > 0]),
        't2w': sum([1 for sess in session_list if len(mybids.get(subject=subject, extension='nii.gz', suffix='T2w', session=sess)) > 0]),
        'flair': sum([1 for sess in session_list if len(mybids.get(subject=subject, extension='nii.gz', suffix='FLAIR', session=sess)) > 0]),
    }
    stats_dict[subject]['total_non_repeated_images'] = sum([v for v in stats_dict[subject]['imaging'].values()])

subject_list_NC = [subject for subject in subject_list if stats_dict[subject]['dx'] == 'NC']
subject_list_PD = [subject for subject in subject_list if stats_dict[subject]['dx'] == 'PD']
print(' * Number of subjects: ' + str(len(subject_list)))
print(' * Number of sessions (total): ' + str(sum([sd['num_sessions'] for s, sd in stats_dict.items()])))
print(' * Number of sessions (per subject): ' + str(np.mean([sd['num_sessions'] for s, sd in stats_dict.items()])) + '+/-' + str(np.std([sd['num_sessions'] for s, sd in stats_dict.items()])))

print(' * Number of subjects by group:', end=' ', flush=True)
print('NC=' + str(len(subject_list_NC)), end=' ', flush=True)
print('PD=' + str(len(subject_list_PD)))
print(' * Number of sessions by group (per subject):', end=' ', flush=True)

print('NC=' + str(np.mean([stats_dict[s]['num_sessions'] for s in subject_list_NC])), end=' +/- ', flush=True)
print('NC=' + str(np.std([stats_dict[s]['num_sessions'] for s in subject_list_NC])), end='\n', flush=True)
print('PD=' + str(np.mean([stats_dict[s]['num_sessions'] for s in subject_list_PD])), end=' +/- ', flush=True)
print('PD=' + str(np.std([stats_dict[s]['num_sessions'] for s in subject_list_PD])))

print(' * Timespan (per subject): ' + str(np.mean([sd['timespan'] for s, sd in stats_dict.items()])) + '+/-' + str(np.std([sd['timespan'] for s, sd in stats_dict.items()])))


import seaborn as sns
import matplotlib.pyplot as plt
from os.path import dirname

dataframe = {
    "age": [],
    "dx": [],
    "timespan": []
}
for sid, s_stats in stats_dict.items():
    dataframe['age'] += [float(s_stats['age_bl'])]
    dataframe['dx'] += [s_stats['dx']]
    dataframe['timespan'] += [float(s_stats['timespan'])]

plt.figure()
sns.histplot(dataframe, x="age", hue="dx", element="step")
plt.savefig(join(dirname(bids_dir), 'age_hist.png'))

plt.figure()
sns.histplot(dataframe, x="timespan", hue="dx", element="step")
plt.savefig(join(dirname(bids_dir), 'timespan_hist.png'))


