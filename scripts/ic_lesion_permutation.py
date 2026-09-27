"""Exploratory exact subtype-label permutation after P7 result; not confirmatory."""
import itertools,sys
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from ubiomark import geo
path=geo.download_matrices('GSE28242')[0];expr,ann,_=geo.parse_series_matrix(path)
state=ann.Sample_characteristics_ch1_2.str.lower()
ctrl=ann.index[state.eq('disease state: normal')]
cases=ann.index[~state.eq('disease state: normal')]
lesion=tuple(ann.index[state.eq('disease state: pbs with lesions')])
x=geo.to_gene_level(expr,geo.probe_to_symbol('GPL6244'))
disc=pd.read_csv('results/meta_discovery/interstitial_cystitis.csv.gz',index_col=0).sort_values('p').head(50)
genes=disc.index.intersection(x.index)
sign=np.sign(disc.loc[genes,'mu'].to_numpy())
X=x.loc[genes]; ref=X.loc[:,ctrl].mean(axis=1).to_numpy()

def stat(lesion_ids):
 a=X.loc[:,list(lesion_ids)].mean(axis=1).to_numpy()-ref
 b=X.loc[:,list(set(cases)-set(lesion_ids))].mean(axis=1).to_numpy()-ref
 return float(np.mean(np.sign(a)==sign)-np.mean(np.sign(b)==sign))
obs=stat(lesion)
allstats=np.array([stat(ids) for ids in itertools.combinations(cases,3)])
assert len(allstats)==56
p=float(np.mean(allstats>=obs-1e-12))
out={'n_genes':len(genes),'observed_concordance_difference':obs,'n_exact_partitions':56,'partitions_at_least_observed':int((allstats>=obs-1e-12).sum()),'exact_one_sided_p':p,'null_mean':float(allstats.mean()),'null_sd':float(allstats.std(ddof=1)),'status':'exploratory post-P7, same samples and shared controls'}
import json
with open('results/ic_lesion_permutation.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
