"""Source-verified, within-cell-line descriptive APF-versus-mock effects.

Not IC-vs-control disease data, eight nested records and only two cell lines.
"""
import csv
import gzip
import hashlib
import io
import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
from ubiomark.geo import parse_series_matrix
PIN='ff6dd2c93d2146406d0021a6a843694be23a5d4a710c82c7d8bb866c9eed6967'
SERIES='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE621'
# Source titles/description independently verified, with per-line APF/mock replication.
GROUPS={'line1_APF':['GSM4877','GSM4879'], 'line1_mock':['GSM4878','GSM4880'],
        'line2_APF':['GSM4881','GSM4883'], 'line2_mock':['GSM4882','GSM4884']}

def run():
    source=ROOT/'data/raw/matrix/GSE621_series_matrix.txt.gz'
    if hashlib.sha256(source.read_bytes()).hexdigest()!=PIN:
        raise ValueError('GSE621 pinned source changed')
    x,ann,_=parse_series_matrix(str(source))
    if x.shape!=(4132,28) or set(ann.Sample_platform_id)!={'GPL262'}:
        raise ValueError('source series shape/platform changed')
    ids=sum(GROUPS.values(),[])
    if len(set(ids))!=8: raise ValueError('sample ID duplication')
    old=list(csv.DictReader(open(ROOT/'results/series/interstitial_cystitis__GSE621.labels.csv')))
    oldlabels={r['gsm']:r['label'] for r in old}
    rows=[]
    for group, gsm_ids in GROUPS.items():
        cell='1' if 'line1' in group else '2'
        treatment='mock' if 'mock' in group else 'APF'
        for gsm in gsm_ids:
            if oldlabels[gsm]: raise ValueError(f'{gsm} unexpectedly in historical IC diagnosis model')
            url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full'
            file=HERE/f'{gsm}.soft.txt.gz'
            if not file.exists():
                with urllib.request.urlopen(url, timeout=35) as r: raw=r.read()
                if len(raw)<35000: raise ValueError('GEO source record too short')
                with gzip.open(file,'wb',compresslevel=9) as f:f.write(raw)
            else: raw=gzip.open(file,'rb').read()
            text=raw.decode().replace('\r\n','\n')
            title=ann.loc[gsm,'Sample_title']
            desc=ann.loc[gsm,'Sample_description_0']
            for key in [f'^SAMPLE = {gsm}',f'!Sample_geo_accession = {gsm}',
                        f'!Sample_title = {title}',f'!Sample_description = {desc}',
                        '!Sample_organism_ch1 = Homo sapiens','!Sample_platform_id = GPL262',
                        '!Sample_series_id = GSE621']:
                if key not in text.splitlines():raise ValueError(f'{gsm}: {key} mismatch')
            if not re.search(rf'cell line\s*{cell}',title+' '+desc,re.I):
                raise ValueError(f'{gsm} cell line mismatch')
            if treatment=='mock' and 'mock' not in title.lower():
                raise ValueError(f'{gsm} missing mock in title')
            if treatment=='APF' and ('mock' in title.lower() or 'apf' not in title.lower()):
                raise ValueError(f'{gsm} treatment mismatch')
            start=text.index('!sample_table_begin\n')+len('!sample_table_begin\n')
            end=text.index('!sample_table_end',start)
            table=pd.read_csv(io.StringIO(text[start:end]),sep='\t',index_col=0)
            table.index=table.index.astype(str)
            if len(table)!=4132 or set(table.index)!=set(x.index.astype(str)):
                raise ValueError(f'{gsm} probe table mismatch')
            a=table.loc[x.index.astype(str),'VALUE'].to_numpy(float)
            b=x[gsm].to_numpy(float)
            if not np.allclose(a,b,atol=1e-7,rtol=1e-7,equal_nan=True):
                raise ValueError(f'{gsm} source/series values differ')
            rows.append({'gsm':gsm,'series':SERIES,'url':url,'sha256':hashlib.sha256(raw).hexdigest(),
                         'title':title,'description':desc,'cell_line':cell,'treatment':treatment,
                         'source_matrix_sha256':PIN,'source_probes_matched':len(table)})
    pd.DataFrame(rows).to_csv(HERE/'GSE621_treatment_crosswalk.csv',index=False)
    # Block on the independent line; two APF and two mock measurements per line.
    d1=x[GROUPS['line1_APF']].mean(axis=1)-x[GROUPS['line1_mock']].mean(axis=1)
    d2=x[GROUPS['line2_APF']].mean(axis=1)-x[GROUPS['line2_mock']].mean(axis=1)
    effects=pd.DataFrame({'probe':x.index.astype(str),'line1_APF_minus_mock':d1.to_numpy(),
                          'line2_APF_minus_mock':d2.to_numpy()})
    effects['mean_APF_minus_mock']=(effects.line1_APF_minus_mock+effects.line2_APF_minus_mock)/2
    effects['same_sign']=(effects.line1_APF_minus_mock*effects.line2_APF_minus_mock)>0
    effects.to_csv(HERE/'GSE621_APF_two_line_probe_effects.csv.gz',index=False,compression='gzip')
    result={'source':SERIES,'matrix_sha256':PIN,'n_used_GSM':len(ids),'n_cell_lines':2,
            'n_APF_libraries':4,'n_mock_libraries':4,'n_probes':len(effects),
            'finite_probe_pairs':int((np.isfinite(d1)&np.isfinite(d2)).sum()),
            'n_effect_sign_agree':int(effects.same_sign.sum()),
            'caveat':'Post-exposure treatment response in two cell lines, not eight independent donors. Cell-line replicates and membrane effects are not resolved. No disease comparison, biomarker or published benchmark win.'}
    (HERE/'GSE621_treatment_result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':print(json.dumps(run(),indent=2))
