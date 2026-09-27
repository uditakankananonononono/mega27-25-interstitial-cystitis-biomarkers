"""Post-audit IC GSE621 effect sensitivity excluding two IC-derived mock-APF controls.
Not a new independent cohort or a replacement of the frozen, previously viewed result.
"""
import sys,json
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from ubiomark import geo,stats
G='GSE621';_,ann,_=geo.parse_series_matrix(geo.download_matrices(G)[0]);gene=pd.read_pickle('data/processed/interstitial_cystitis__GSE621.pkl.gz');old=pd.read_csv(f'results/series/interstitial_cystitis__{G}.labels.csv').set_index('gsm').label
case=ann.Sample_title.str.fullmatch(r'ic patient ?[1-6]',case=False)
ctrl=ann.Sample_title.str.fullmatch(r'normal control [1-6]',case=False)
assert (case.sum(),ctrl.sum())==(6,6)
assert set(old[old=='case'].index)==set(ann.index[case])
assert set(old[old=='control'].index)-set(ann.index[ctrl])=={'GSM4869','GSM4872'}
assert all('mock' in ann.loc[g,'Sample_title'].lower() and 'ic' in ann.loc[g,'Sample_title'].lower() for g in ('GSM4869','GSM4872'))
assert set(gene.columns)==set(old[old.notna()].index) # use previously verified matrix GSMs; no reparsing drift
x1=gene.loc[:,ann.index[case]].to_numpy(float);x0=gene.loc[:,ann.index[ctrl]].to_numpy(float)
g,v=stats.hedges_g(x1,x0)
out=pd.DataFrame({'gene':gene.index,'g':g,'v':v}).replace([np.inf,-np.inf],np.nan).dropna()
out.to_csv('results/ic_GSE621_strict_controls_effects.csv.gz',index=False)
oldfx=pd.read_csv('results/series/interstitial_cystitis__GSE621.csv.gz').set_index('gene');new=out.set_index('gene');shared=new.index.intersection(oldfx.index)
diff=(new.loc[shared,'g']-oldfx.loc[shared,'g']).abs()
res={'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE621','type':'post-audit effect-only phenotype correction; not a new held-out test','case_ic_patients':6,'normal_controls':6,'old_controls':8,'excluded_old_controls':['GSM4869','GSM4872'],'old_controls_problem':'IC-derived mock-APF sample titles, not normal donors','genes_old':len(oldfx),'genes_strict':len(new),'genes_shared':len(shared),'median_absolute_hedges_g_shift':float(diff.median()),'n_sign_changed':int((np.sign(new.loc[shared,'g'])!=np.sign(oldfx.loc[shared,'g'])).sum())}
with open('results/ic_GSE621_strict_controls_sensitivity.json','w') as f:json.dump(res,f,indent=2)
print(json.dumps(res,indent=2))
