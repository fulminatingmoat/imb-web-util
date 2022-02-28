import gc
import hashlib
import os
import shutil

import celery.exceptions
import numpy as np
import pandas as pd
import scipy.stats as st

COMPRESSED_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'compressed')
TASK_RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'task_results')
RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'results')

MATERNAL_REQUIRED_HEADERS = ['rsid', 'chr', 'pos', 'ea', 'nea', 'eaf', 'beta', 'se', 'p']
FETAL_REQUIRED_HEADERS = ['rsid', 'chr', 'pos', 'ea', 'nea', 'eaf', 'beta', 'se', 'p']

default_maternal = '7b59a43b06a6a3d507321dd05a4b5bc4'
default_fetal = '10fe5ad6852fd0fb1f40c77717f33151'


def optimize_floats(df: pd.DataFrame) -> pd.DataFrame:
    floats = df.select_dtypes(include=['float64']).columns.tolist()
    df[floats] = df[floats].apply(pd.to_numeric, downcast='float')
    return df


def optimize_ints(df: pd.DataFrame) -> pd.DataFrame:
    ints = df.select_dtypes(include=['int64']).columns.tolist()
    df[ints] = df[ints].apply(pd.to_numeric, downcast='integer')
    return df


def optimize_objects(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.select_dtypes(include=['object']):
        if not (type(df[col][0]) == list):
            num_unique_values = len(df[col].unique())
            num_total_values = len(df[col])
            if float(num_unique_values) / num_total_values < 0.5:
                df[col] = df[col].astype('category')
    return df


def optimize(df: pd.DataFrame):
    return optimize_floats(optimize_ints(optimize_objects(df)))


def impute(self, maternal_file=None, fetal_file=None, intercept='0.1352'):
    # Read in all the data
    self.update_state(state='PROGRESS', meta={'progress': 0.1, 'total': 1, 'status': 'Loading maternal file'})
    maternal_file = pd.read_pickle(os.path.join(COMPRESSED_DIR, f"{maternal_file}.pkl"))
    maternal_file.drop(columns=[col for col in maternal_file if col not in MATERNAL_REQUIRED_HEADERS], inplace=True)
    self.update_state(state='PROGRESS', meta={'progress': 0.15, 'total': 1, 'status': 'Loading fetal file'})
    fetal_file = pd.read_pickle(os.path.join(COMPRESSED_DIR, f"{fetal_file}.pkl"))
    fetal_file.drop(columns=[col for col in fetal_file if col not in FETAL_REQUIRED_HEADERS], inplace=True)

    # Check that the headers are correct

    for x in MATERNAL_REQUIRED_HEADERS:
        if x not in maternal_file.columns:
            self.update_state(state='FAILURE',
                              meta={'progress': 0.5, 'total': 1,
                                    'status': f'Maternal file is missing the {x} header {maternal_file.columns}',
                                    'exc_type': 'KeyError',
                                    'exc_message': f'Maternal file is missing the {x} header {maternal_file.columns}'})
            raise celery.exceptions.Ignore()

    for x in FETAL_REQUIRED_HEADERS:
        if x not in fetal_file.columns:
            self.update_state(state='FAILURE',
                              meta={'progress': 0.5, 'total': 1,
                                    'status': f'Fetal file is missing the {x} header {fetal_file.columns}',
                                    'exc_type': 'KeyError',
                                    'exc_message': f'Fetal file is missing the {x} header {fetal_file.columns}'})
            raise celery.exceptions.Ignore()

    intercept = float(intercept)

    maternal_file.rename(
        columns={'ea': 'ea_maternal', 'nea': 'nea_maternal', 'eaf': 'eaf_maternal', 'beta': 'beta_maternal',
                 'se': 'se_maternal', 'p': 'p_maternal'}, inplace=True)
    fetal_file.rename(
        columns={'ea': 'ea_fetal', 'nea': 'nea_fetal', 'eaf': 'eaf_fetal', 'beta': 'beta_fetal', 'se': 'se_fetal',
                 'p': 'p_fetal'}, inplace=True)
    self.update_state(state='PROGRESS', meta={'progress': 0.20, 'total': 1, 'status': 'Merging data'})

    merged_data = pd.merge(maternal_file, fetal_file, on=['rsid', 'chr', 'pos'])

    print(merged_data.info())

    del maternal_file
    del fetal_file
    print(gc.collect())

    optimize(merged_data)
    print(merged_data.info())

    self.update_state(state='PROGRESS', meta={'progress': 0.25, 'total': 1, 'status': 'Matching columns'})

    merged_data['beta_maternal'] = merged_data.apply(lambda a: a['beta_maternal'] * -1 if a['ea_fetal'] == a[
        'nea_maternal'] else a['beta_maternal'], axis=1)
    merged_data['ea_maternal'] = merged_data.apply(lambda a: a['nea_maternal'] if a['ea_fetal'] == a[
        'nea_maternal'] else a['ea_maternal'], axis=1)
    merged_data['nea_maternal'] = merged_data.apply(lambda a: a['ea_maternal'] if a['ea_fetal'] == a[
        'nea_maternal'] else a['nea_maternal'], axis=1)
    merged_data['eaf_maternal'] = merged_data.apply(lambda a: 1 - a['eaf_maternal'] if a['ea_fetal'] == a[
        'nea_maternal'] else a['eaf_maternal'], axis=1)

    optimize(merged_data)

    self.update_state(state='PROGRESS', meta={'progress': 0.30, 'total': 1, 'status': 'Adjusting'})

    merged_data['fetal_beta_adjusted'] = ((4 / 3) * merged_data['beta_fetal']) - (
            (2 / 3) * merged_data['beta_maternal'])
    merged_data['maternal_beta_adjusted'] = ((4 / 3) * merged_data['beta_maternal']) - (
            (2 / 3) * merged_data['beta_fetal'])
    merged_data['fetal_var_adjusted'] = (16 / 9) * merged_data['se_fetal'] ** 2 + (4 / 9) * merged_data[
        'se_maternal'] ** 2 - (16 / 9) * merged_data['se_maternal'] * merged_data['se_fetal'] * intercept
    merged_data['maternal_var_adjusted'] = (16 / 9) * merged_data['se_maternal'] ** 2 + (4 / 9) * merged_data[
        'se_fetal'] ** 2 - (16 / 9) * merged_data['se_maternal'] * merged_data['se_fetal'] * intercept

    merged_data['fetal_se_adjusted'] = merged_data['fetal_var_adjusted'] ** 0.5
    merged_data['maternal_se_adjusted'] = merged_data['maternal_var_adjusted'] ** 0.5

    merged_data['covar'] = (20 / 9) * merged_data['se_maternal'] * merged_data['se_fetal'] * intercept - (8 / 9) * \
                           merged_data['se_maternal'] ** 2 - (8 / 9) * merged_data['se_fetal'] ** 2

    optimize(merged_data)
    self.update_state(state='PROGRESS', meta={'progress': 0.35, 'total': 1, 'status': 'ChiSq'})
    print(merged_data.columns)
    print(merged_data.info())
    gc.collect()
    merged_data['chisq2df'] = merged_data.apply(lambda a: (np.array(
        [a['fetal_beta_adjusted'], a['maternal_beta_adjusted']]) @ np.linalg.inv(np.array(
        [[a['fetal_var_adjusted'], a['covar']],
         [a['covar'], a['maternal_var_adjusted']]])) @ np.transpose(
        np.array([[a['fetal_beta_adjusted'], a['maternal_beta_adjusted']]])))[0], axis=1)

    print(merged_data.head())

    self.update_state(state='PROGRESS', meta={'progress': 0.35, 'total': 1, 'status': 'ChiSq adjusting'})
    print('adjusting')
    merged_data['pval2df'] = st.chi2.sf(merged_data['chisq2df'], 2)
    merged_data["chisq_fetal_adj"] = merged_data.apply(
        lambda a: a['fetal_beta_adjusted'] ** 2 / a['fetal_var_adjusted'], axis=1)
    merged_data["chisq_maternal_adj"] = merged_data.apply(
        lambda a: a['maternal_beta_adjusted'] ** 2 / a['maternal_var_adjusted'], axis=1)
    merged_data["pval_fetal_adj"] = st.chi2.sf(merged_data['chisq_fetal_adj'], 1)
    merged_data['pval_maternal_adj'] = st.chi2.sf(merged_data['chisq_maternal_adj'], 1)
    optimize(merged_data)
    print('adjusted')
    results = {'CHR': merged_data['chr'], 'SNP': merged_data['rsid'], 'BP': merged_data['pos'],
               'A1': merged_data['ea_fetal'], 'A2': merged_data['nea_fetal'], 'FREQ': merged_data['eaf_fetal'],
               'BETA_F': merged_data['beta_fetal'], 'SE_F': merged_data['se_fetal'], "PVAL_F": merged_data['p_fetal'],
               "BETA_M": merged_data['beta_maternal'], "SE_M": merged_data['se_maternal'],
               "PVAL_M": merged_data['p_maternal'], "BETA_F_ADJ": merged_data['fetal_beta_adjusted'],
               "PVAL_F_ADJ": merged_data['pval_fetal_adj'], "BETA_M_ADJ": merged_data['maternal_beta_adjusted'],
               "PVAL_M_ADJ": merged_data['pval_maternal_adj'],
               "CHISQ_2DF": merged_data['chisq2df'], "PVAL_2DF": merged_data['pval2df'],
               "F_SE_ADJ": merged_data['fetal_se_adjusted'],
               "M_SE_ADJ": merged_data['maternal_se_adjusted']}

    del merged_data
    gc.collect()

    self.update_state(state='PROGRESS', meta={'progress': 0.6, 'total': 1, 'status': 'Imputing'})

    df = pd.DataFrame(results)

    filename = f"{self.request.id}.tsv"

    df.to_csv(os.path.join(TASK_RESULT_DIR, filename), '\t', index=False)

    md5_hash = hashlib.md5()
    with open(os.path.join(TASK_RESULT_DIR, filename), "rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)

    if not os.path.exists(os.path.join(RESULT_DIR, f"{md5_hash.hexdigest()}.tsv")):
        shutil.move(os.path.join(TASK_RESULT_DIR, filename), os.path.join(RESULT_DIR, f"{md5_hash.hexdigest()}.tsv"))
    else:
        os.remove(os.path.join(TASK_RESULT_DIR, filename))

    self.update_state(state='SUCCESS',
                      meta={'progress': 1, 'total': 1, 'status': 'Done', 'result': md5_hash.hexdigest()})

    return {'filename': f"{filename}.tsv", 'md5_hash': md5_hash.hexdigest(), 'progress': 1, 'total': 1,
            'status': 'Done', 'result': md5_hash.hexdigest()}
