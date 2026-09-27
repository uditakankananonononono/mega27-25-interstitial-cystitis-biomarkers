"""Post-audit GSE11783 donor-level IC sensitivity: average paired biopsies within five IC donors.
Historical sample-level estimates remain available, not silently overwritten.
"""
import re,sys,json,collections
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from ubiomark import geo,stats
G='GSE11783';_,ann,_=geo.parse_series_matrix(geo.download_matrices(G)[0]);x=pd.read_pickle(f'data/processed/interstitial_cystitis__{G}.pkl.gz');lab=pd.read_csv(f'results/series/interstitial_cystitis__{G}.labels.csv').set_index('gsm').label
assert set(x.columns)==set(lab.index) and len(x.columns)==16
person={g:re.search(r'patient:(\d+)',ann.loc[g,'Sample_title']).group(1) for g in x.columns}
group={p:lab.loc[[g for g in x.columns if person[g]==p]].unique().tolist() for p in set(person.values())};assert all(len(v)==1 for v in group.values())
assert collections.Counter(v[0] for v in group.values())=={'case':5,'control':6}
pooled=pd.DataFrame({p:x[[g for g in x.columns if person[g]==p]].mean(axis=1) for p in sorted(group)})
cases=[p for p,v in group.items() if v==['case']];ctrl=[p for p,v in group.items() if v==['control']]
g,v=stats.hedges_g(pooled[cases].to_numpy(float),pooled[ctrl].to_numpy(float))
out=pd.DataFrame({'gene':pooled.index,'g':g,'v':v}).replace([np.inf,-np.inf],np.nan).dropna();out.to_csv('results/ic_GSE11783_donor_effects.csv.gz',index=False)
old=pd.read_csv('results/series/interstitial_cystitis__GSE11783.csv.gz').set_index('gene');new=out.set_index('gene');ix=old.index.intersection(new.index)
result={'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11783','type':'post-audit donor-level effect sensitivity, not new validation','n_case_libraries':10,'n_case_people':5,'n_control_libraries':6,'n_control_people':6,'case_people_two_biopsies':5,'genes_shared':len(ix),'median_absolute_hedges_g_shift':float((new.loc[ix,'g']-old.loc[ix,'g']).abs().median()),'n_effect_sign_changed':int((np.sign(new.loc[ix,'g'])!=np.sign(old.loc[ix,'g'])).sum()),'control':'historical sample-level effects retained; no downstream discovery/benchmark recomputation'}
with open('results/ic_GSE11783_donor_sensitivity.json','w') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
