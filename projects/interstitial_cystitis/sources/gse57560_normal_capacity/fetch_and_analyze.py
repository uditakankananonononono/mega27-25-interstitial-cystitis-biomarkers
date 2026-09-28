"""Fetch and verify previously omitted GSE57560 normal-capacity IC biopsy records.

Exploratory capacity-stratified probe analysis, same series as old tiny holdout;
not independent validation, not a clinical marker or corrected gene-level meta.
"""
import csv
import gzip
import hashlib
import io
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from ubiomark.geo import parse_series_matrix
from ubiomark.stats import hedges_g, bh_fdr

URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57560/matrix/GSE57560_series_matrix.txt.gz'
PIN = 'bdd66a6e3be02701af2fcbc06695dead3f231845938f194a1d04a5dfaf71f255'


def run():
    matrix_file = ROOT / 'data/geo/gse57560/GSE57560_series_matrix.txt.gz'
    if hashlib.sha256(matrix_file.read_bytes()).hexdigest() != PIN:
        raise ValueError('GSE57560 GEO matrix digest changed')
    x, ann, _ = parse_series_matrix(str(matrix_file))
    if x.shape != (62976, 16) or set(ann.Sample_platform_id) != {'GPL16699'}:
        raise ValueError('unexpected GEO matrix dimensions/platform')
    original = list(csv.DictReader(open(ROOT / 'results/series/interstitial_cystitis__GSE57560.labels.csv')))
    historical = {r['gsm']:r['label'] for r in original if r['label']}
    if len(historical) != 7 or any(g not in ann.index for g in historical):
        raise ValueError('historical seven-library source mismatch')
    rows=[]; groups={'normal_capacity_IC':[], 'low_capacity_IC':[], 'healthy_control':[]}
    for gsm, a in ann.iterrows():
        char = [a[f'Sample_characteristics_ch1_{i}'] for i in range(4)]
        if 'phenotype: normal capactiy' in char and 'disease state: Interstitial cystitis' in char:
            group = 'normal_capacity_IC'
        elif 'phenotype: low capacity' in char and 'disease state: Interstitial cystitis' in char:
            group = 'low_capacity_IC'
        elif 'phenotype: control' in char and 'disease state: control' in char:
            group = 'healthy_control'
        else:
            raise ValueError(f'unknown clinical stratum for {gsm}: {char}')
        expected = {'normal_capacity_IC':'', 'low_capacity_IC':'case', 'healthy_control':'control'}[group]
        if next(r['label'] for r in original if r['gsm']==gsm) != expected:
            raise ValueError(f'old label mismatch for {gsm}')
        groups[group].append(gsm)
        # Individually source-check the nine previously excluded specimens.
        if group == 'normal_capacity_IC':
            url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full'
            file=HERE/f'{gsm}.soft.txt.gz'
            if not file.exists():
                with urllib.request.urlopen(url, timeout=35) as response:
                    raw=response.read()
                if len(raw)<10000:
                    raise ValueError(f'source record too short {gsm}')
                with gzip.open(file,'wb',compresslevel=9) as out: out.write(raw)
            else:
                raw=gzip.open(file,'rb').read()
            text=raw.decode()
            for line in [f'^SAMPLE = {gsm}', f'!Sample_geo_accession = {gsm}',
                         f'!Sample_title = {a.Sample_title}',
                         '!Sample_platform_id = GPL16699',
                         '!Sample_organism_ch1 = Homo sapiens',
                         '!Sample_characteristics_ch1 = phenotype: normal capactiy',
                         '!Sample_characteristics_ch1 = disease state: Interstitial cystitis']:
                if line not in text.splitlines():
                    raise ValueError(f'{gsm} source field mismatch {line}')
            text=text.replace("\r\n","\n")
            start=text.index('!sample_table_begin\n')+len('!sample_table_begin\n')
            end=text.index('!sample_table_end',start)
            t=pd.read_csv(io.StringIO(text[start:end]),sep='\t',index_col=0)
            t.index=t.index.astype(str)
            if not set(x.index.astype(str)) == set(t.index.astype(str)) or len(t)!=62976:
                raise ValueError(f'{gsm} probe set mismatch')
            y=t.loc[x.index.astype(str),'VALUE'].to_numpy(dtype=float)
            z=x[gsm].to_numpy(dtype=float)
            # float32 conversion in parser imposes a small relative error.
            if not np.allclose(y,z,atol=.0002,rtol=1e-7,equal_nan=True):
                raise ValueError(f'{gsm} source probe values differ from series matrix')
            rows.append({'gsm':gsm,'group':group,'source_url':url,'source_sha256':hashlib.sha256(raw).hexdigest(),
                         'title':a.Sample_title,'source_phenotype':char[2],'source_disease':char[3],
                         'matrix_sha256':PIN,'probes_matched':len(t)})
    if {k:len(v) for k,v in groups.items()} != {'normal_capacity_IC':9,'low_capacity_IC':4,'healthy_control':3}:
        raise ValueError('unexpected clinical group counts')
    pd.DataFrame(rows).to_csv(HERE/'GSE57560_normal_capacity_used_crosswalk.csv',index=False)
    # Contrast is disease severity/capacity within IC; never label normal-capacity IC as healthy controls.
    g, var = hedges_g(x[groups['low_capacity_IC']].to_numpy(dtype=float),
                      x[groups['normal_capacity_IC']].to_numpy(dtype=float))
    _, pvals=ttest_ind(x[groups['low_capacity_IC']].to_numpy(dtype=float),
                        x[groups['normal_capacity_IC']].to_numpy(dtype=float),
                        axis=1, equal_var=False, nan_policy='omit')
    q=bh_fdr(pvals)
    out=pd.DataFrame({'probe':x.index.astype(str),'g_low_minus_normal':g,'v':var,'welch_p':pvals,'BH_q':q})
    out.to_csv(HERE/'GSE57560_capacity_probe_effects.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    result={'source':URL,'source_sha256':PIN,'n_normal_capacity_IC':9,'n_low_capacity_IC':4,
            'n_healthy_controls_excluded':3,'n_source_records_newly_analyzed':len(rows),
            'n_probes':len(out),'n_finite_effects':int(np.isfinite(g).sum()),
            'welch_BH_q_below_05':int(np.sum(q<.05)),
            'caveat':'Post-exposure same-cohort capacity-stratified probe sensitivity. Historical 4-vs-3 low-capacity IC versus healthy-control holdout is already seen. These nine are IC patients, not healthy controls. Exploratory per-probe Welch tests without a prespecified endpoint or batch adjustment are not clinical validation; no independent sample, clinical classifier, gene-level panel or published benchmark win.'}
    (HERE/'GSE57560_capacity_result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__': print(json.dumps(run(),indent=2))
