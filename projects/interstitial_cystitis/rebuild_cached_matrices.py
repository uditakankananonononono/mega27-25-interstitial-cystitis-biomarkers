"""Rebuild two exposed IC gene matrices from pinned GEO series/platform sources."""
import gzip, hashlib, io, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from ubiomark.geo import parse_series_matrix,_clean_symbol,to_gene_level

SOURCES={
 'GSE621':{'gpl':'GPL262','matrix_sha256':'ff6dd2c93d2146406d0021a6a843694be23a5d4a710c82c7d8bb866c9eed6967','platform_sha256':'54a078fb563741608f05c347ce3612b34ff7df1616648e8dda25d6b39abcf669','cached_sha256':'3a2849b561ceda924c2de8423753c33fd67e806db970b8507e3d72393147f01d'},
 'GSE11783':{'gpl':'GPL570','matrix_sha256':'05f4b1f257cb5ce5d63584930ba215828356d73bf79d35bbdc38f3c528beae12','platform_sha256':'d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394','cached_sha256':'015d08d24fbce47304a5bbf0023e5997f7370fc686b82ce167acc4384ff189e6'}}

def rebuild(acc):
 c=SOURCES[acc];raw=ROOT/f'data/raw/matrix/{acc}_series_matrix.txt.gz';gpl=ROOT/f'data/raw/gpl/{c["gpl"]}.annot.gz';cache=ROOT/f'data/processed/interstitial_cystitis__{acc}.pkl.gz'
 for name,path in [('matrix',raw),('platform',gpl),('cached',cache)]:
  if hashlib.sha256(path.read_bytes()).hexdigest()!=c[name+'_sha256']:
   raise ValueError(f'{name} SHA-256 differs from pinned source')
 expr,ann,_=parse_series_matrix(str(raw))
 if set(ann.Sample_platform_id)!={c['gpl']}:
  raise ValueError('unexpected platform')
 with gzip.open(gpl,'rt',errors='replace') as f:rows=[line for line in f if not line.startswith(('^','!','#'))]
 t=pd.read_csv(io.StringIO(''.join(rows)),sep='\t',dtype=str,low_memory=False)
 p2s=t.set_index('ID')['Gene symbol'].map(lambda s:_clean_symbol(s,'Gene symbol')).dropna();p2s.index=p2s.index.astype(str)
 gene=to_gene_level(expr,p2s);old=pd.read_pickle(cache)
 if set(gene.index)!=set(old.index) or not set(old.columns)<=set(gene.columns):
  raise ValueError('cached genes or GSMs differ from reconstructed source')
 a=gene.loc[old.index,old.columns].to_numpy();b=old.to_numpy()
 delta=np.abs(a-b)
 threshold=0 if acc=='GSE621' else 1e-6
 if not np.allclose(a,b,atol=threshold,rtol=0,equal_nan=True):
  raise ValueError(f'cached value mismatch > {threshold}')
 return {'study':acc,'raw_probe_by_GSM':list(expr.shape),'cached_gene_by_GSM':list(old.shape),
         'source_url':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}',
         'max_absolute_value_delta':float(np.nanmax(delta)),'absolute_tolerance':threshold,
         'gene_id_set_and_selected_GSMs_reproduced':True,'source_hashes':c}

if __name__=='__main__':
 r=[rebuild(x) for x in SOURCES]
 out=ROOT/'projects/interstitial_cystitis/source_matrix_rebuild_audit.json';out.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps(r,indent=2))
