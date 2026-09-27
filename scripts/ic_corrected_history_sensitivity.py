"""Post-audit IC sensitivity replacing contaminated GSE621 controls and GSE11783 sample-level variance.
Not a fresh validation: GSE57560 labels and outcomes had already been inspected; retain frozen outputs.
"""
import json,sys
import numpy as np,pandas as pd
from scipy.stats import norm
sys.path.insert(0,'src')
from ubiomark import stats
A=pd.read_csv('results/ic_GSE621_strict_controls_effects.csv.gz').set_index('gene').add_suffix('_a')
B=pd.read_csv('results/ic_GSE11783_donor_effects.csv.gz').set_index('gene').add_suffix('_b')
F=pd.concat([A,B],axis=1)
r=stats.dersimonian_laird(F[['g_a','g_b']].to_numpy(float),F[['v_a','v_b']].to_numpy(float))
meta=pd.DataFrame({k:r[k] for k in ['mu','se','z','p','tau2','Q','I2','k']},index=F.index)
meta=meta[meta.k>=2];meta['q']=stats.bh_fdr(meta.p.to_numpy());meta.index.name='gene';meta.sort_values('p').to_csv('results/ic_post_audit_discovery_meta_sensitivity.csv.gz')
old=pd.read_csv('results/meta_discovery/interstitial_cystitis.csv.gz',index_col=0)
val=pd.read_csv('results/series/interstitial_cystitis__GSE57560.csv.gz').set_index('gene')
val['z']=val.g/np.sqrt(val.v);val['p']=2*norm.sf(np.abs(val.z))
common=meta.index.intersection(val.index);top=meta.sort_values('p').head(50);sel=top.index.intersection(common)
s=np.sign(meta.loc[sel,'mu'].to_numpy());v=val.loc[sel];agree=(np.sign(v.z.to_numpy())==s);hits=int((agree&(v.p.to_numpy()<.05)).sum())
# Reuse historical random-gene sign null method, but explicitly post-audit because top list changed.
rng=np.random.default_rng(20260925);Bnull=10000;n=len(sel);cg=np.array(common);ds=np.sign(meta.loc[cg,'mu'].to_numpy());z=val.loc[cg,'z'].to_numpy();pv=val.loc[cg,'p'].to_numpy();n1=np.empty(Bnull);n2=np.empty(Bnull)
for i in range(Bnull):
 ix=rng.choice(len(cg),n,replace=False);a=np.sign(z[ix])==ds[ix];n1[i]=a.mean();n2[i]=int(np.sum(a&(pv[ix]<.05)))
original_top=old.sort_values('p').head(50).index
res={'audit_type':'post-audit sensitivity using previously seen validation; cannot claim fresh independent replication','validation_source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57560','discovery_sources':['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE621','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11783'],'n_discovery_genes_present_in_both':len(meta),'old_top50_intersection_new_top50':len(set(original_top)&set(top.index)),'n_new_top50_measured':n,'n_val_case_libraries':4,'n_val_control_libraries':3,'n_sign_agree':int(agree.sum()),'fraction_sign_agree':float(agree.mean()),'null_sign_mean':float(n1.mean()),'empirical_p_sign':float((1+(n1>=agree.mean()).sum())/(Bnull+1)),'n_direction_nominal_p05':hits,'null_hit_mean':float(n2.mean()),'empirical_p_hits':float((1+(n2>=hits).sum())/(Bnull+1)),'old_frozen_T1':'40/40 signs in prior 6-vs-8 and 10-vs-6 sample-level discovery; old result left intact. This sensitivity requires evidence from BOTH corrected discovery studies (k=2).' ,'caveat':'A 2025 publication used overlapping studies; a 4-vs-3 validation is too small to establish novelty or clinical classifier utility.'}
with open('results/ic_post_audit_replication_sensitivity.json','w') as f:json.dump(res,f,indent=2)
print(json.dumps(res,indent=2))
