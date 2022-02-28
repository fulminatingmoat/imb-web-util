import hashlib
import os
import shutil

import celery.exceptions
import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.stats as st

COMPRESSED_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'compressed')
TASK_RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'task_results')
RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'results')


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


class LDMapReferencePanel:
    def __init__(self, ld_df, map_df):
        self.ld_df = ld_df
        self.map_df = map_df


DEFAULT_MAF = 0.005
REFERENCE_MAP_REQUIRED_HEADERS = ['SNP_pos', 'SNP_id', "Major", "Minor", 'maf']
SUMMARY_REQUIRED_HEADERS = ['SNP_id', 'SNP_pos', 'Effect_allele', 'Non_effect_allele', 'Z']


def impute(self, maf=DEFAULT_MAF, reference_file=None, map_file=None, summary_file=None,
           genome_build='POS18', custom=False):
    # Read in all the data
    self.update_state(state='PROGRESS', meta={'progress': 0.1, 'total': 1, 'status': 'Loading reference panel'})
    reference_file = pd.read_pickle(os.path.join(COMPRESSED_DIR, f"{reference_file}.pkl"))
    optimize(reference_file)
    self.update_state(state='PROGRESS', meta={'progress': 0.15, 'total': 1, 'status': 'Loading reference map'})
    map_file = pd.read_pickle(os.path.join(COMPRESSED_DIR, f"{map_file}.pkl"))
    optimize(map_file)
    self.update_state(state='PROGRESS', meta={'progress': 0.2, 'total': 1, 'status': 'Loading summary stats'})
    summary_file = pd.read_pickle(os.path.join(COMPRESSED_DIR, f"{summary_file}.pkl"))
    optimize(summary_file)
    if not custom:
        map_file.rename(columns={genome_build: 'SNP_pos'}, inplace=True)
    summary = summary_file

    reference_panel = LDMapReferencePanel(reference_file, map_file)

    for x in REFERENCE_MAP_REQUIRED_HEADERS:
        if x not in map_file.columns:
            self.update_state(state='FAILURE',
                              meta={'progress': 0.5, 'total': 1, 'status': f'Reference map is missing the {x} header',
                                    'exc_type': 'KeyError',
                                    'exc_message': f'Reference map is missing the {x} header'})
            raise celery.exceptions.Ignore()

    for x in SUMMARY_REQUIRED_HEADERS:
        if x not in summary.columns:
            self.update_state(state='FAILURE',
                              meta={'progress': 0.5, 'total': 1, 'status': f'Summary stats is missing the {x} header',
                                    'exc_type': 'KeyError',
                                    'exc_message': f'Summary stats is missing the {x} header'})
            raise celery.exceptions.Ignore()

    self.update_state(state='PROGRESS', meta={'progress': 0.25, 'total': 1, 'status': 'Normalising summary'})

    # Filter out SNPs with MAF < 0.005, or given MAF
    map_maf_filtered = reference_panel.map_df.query('@maf <= maf <= 1-@maf')
    ld_maf_filtered = reference_panel.ld_df.iloc[:, map_maf_filtered.index].loc[map_maf_filtered.index]
    ld_maf_filtered_reset_index = ld_maf_filtered.reset_index()

    # Merge the summary and reference panel
    reference_summary_merged = map_maf_filtered.merge(summary, on="SNP_pos")
    # Check that the alleles match
    reference_summary_matched_positions = reference_summary_merged.query(
        '(Major == Effect_allele & Minor == Non_effect_allele) | (Major == Non_effect_allele & Minor == Effect_allele)')[
        'SNP_pos']

    map_maf_filtered_reset_index = map_maf_filtered.reset_index()
    # get array of positions in summary
    typed_idx = map_maf_filtered_reset_index.loc[
        map_maf_filtered_reset_index['SNP_pos'].isin(reference_summary_matched_positions)].index
    # Use reference to create typed_typed, and typed_untyped matrices
    if len(typed_idx) >= 1:
        ld_tt = ld_maf_filtered.iloc[:, typed_idx]
        ld_tt = ld_tt.iloc[typed_idx, :]
        ld_tt = ld_tt + np.identity(len(ld_tt.columns)) * 0.15
        ld_it = ld_maf_filtered.drop(typed_idx, axis=0, errors='ignore')
        ld_it = ld_it.iloc[:, typed_idx]
        map_t = map_maf_filtered.iloc[typed_idx]
        map_i = map_maf_filtered.drop(typed_idx, axis=0, errors='ignore')
    else:
        self.update_state(state='FAILURE',
                          meta={'progress': 0.5, 'total': 1,
                                'status': 'No matched SNP. Please check your summary file\'s SNP position',
                                'exc_type': 'KeyError',
                                'exc_message': f'No matched SNP. Please check your summary file\'s SNP position'})
        raise celery.exceptions.Ignore()

    typed_in_both_idx = summary.loc[summary['SNP_pos'].isin(reference_summary_matched_positions)].index

    if len(typed_in_both_idx) >= 1:
        summary = summary.iloc[typed_in_both_idx]
    else:
        self.update_state(state='FAILURE',
                          meta={'progress': 0.5, 'total': 1, 'status': 'Could not find SNP positions in summary',
                                'exc_type': 'KeyError',
                                'exc_message': f'Could not find SNP positions in summary'})
        raise celery.exceptions.Ignore()

    Zt = np.full(len(summary['Z']), None)

    self.update_state(state='PROGRESS', meta={'progress': 0.3, 'total': 1, 'status': 'Imputing'})

    for i in range(len(Zt)):
        if summary.iloc[i]['Effect_allele'] != \
                map_t.loc[map_t[map_t['SNP_pos'] == summary.iloc[i]['SNP_pos']].index].iloc[0]['Minor']:
            Zt[map_t['SNP_pos'] == summary.iloc[i]['SNP_pos']] = -summary.iloc[i]['Z']
        else:
            Zt[map_t['SNP_pos'] == summary.iloc[i]['SNP_pos']] = summary.iloc[i]['Z']

    ld_tt_inv = la.solve(ld_tt, np.identity(len(ld_tt.columns)))

    weight = np.dot(ld_it, ld_tt_inv)

    Zi = np.dot(weight, Zt)

    var = np.full(len(weight), None)

    self.update_state(state='PROGRESS', meta={'progress': 0.4, 'total': 1, 'status': 'Imputing'})

    for i in range(len(weight)):
        var[i] = np.sum(np.multiply(np.outer(weight[i], weight[i]), ld_tt.to_numpy()))

    impute_p = lambda x: 2 * st.norm.cdf(-abs(x))
    imputed_p = impute_p(Zi.astype('float64'))

    self.update_state(state='PROGRESS', meta={'progress': 0.5, 'total': 1, 'status': 'Imputing'})

    results = {'Marker_id': map_i['SNP_id'], 'Marker_pos': map_i['SNP_pos'], 'Effect_allele': map_i['Minor'],
               'Non_effect_allele': map_i['Major'], 'Imputed_Z': Zi, 'r2pred': var,
               'imputed_P': 2 * st.norm.cdf(-abs(Zi.astype('float64')))}

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
