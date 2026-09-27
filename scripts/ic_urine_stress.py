"""P7 registered cross-tissue IC/BPS sign-set stress test on GSE28242."""
import os, sys, json, csv
import numpy as np
import pandas as pd
sys.path.insert(0,'src')
from ubiomark import geo, stats

GSE='GSE28242'; URL=geo.series_dir(GSE)+'GSE28242_series_matrix.txt.gz'
path=geo.download_matrices(GSE)[0]
expr,ann,_=geo.parse_series_matrix(path)
assert len(ann)==13 and expr.shape[1]==13 and ann.index.is_unique
assert set(ann.Sample_platform_id)=={'GPL6244'}
state=ann.Sample_characteristics_ch1_2.str.lower()
normal=state.eq('disease state: normal')
without=state.eq('disease state: pbs without lesions')
withlesion=state.eq('disease state: pbs with lesions')
assert (normal.sum(),without.sum(),withlesion.sum())==(5,5,3)
assert all('control' in v for v in ann.loc[normal,'Sample_title'].str.lower())
assert all('without lesions' in v for v in ann.loc[without,'Sample_title'].str.lower())
assert all('with lesions' in v for v in ann.loc[withlesion,'Sample_title'].str.lower())
assert not (set(ann.index) & set(pd.concat([pd.read_csv(f'results/series/interstitial_cystitis__{t}.labels.csv') for t in ('GSE621','GSE11783','GSE57560')]).gsm))
map_=geo.probe_to_symbol('GPL6244')
gene=geo.to_gene_level(expr,map_)
disc=pd.read_csv('results/meta_discovery/interstitial_cystitis.csv.gz',index_col=0)
selected=disc.sort_values('p').head(50)
rows=[]; gene_rows=[]
for stratum,cases in [('all_PBS',without|withlesion),('without_lesions',without),('with_lesions',withlesion)]:
    X1=gene.loc[:,ann.index[cases]].to_numpy(float);X0=gene.loc[:,ann.index[normal]].to_numpy(float)
    g,v=stats.hedges_g(X1,X0)
    effects=pd.DataFrame({'g':g,'v':v},index=gene.index).replace([np.inf,-np.inf],np.nan).dropna()
    common=disc.index.intersection(effects.index)
    observed=selected.index.intersection(common)
    signs=np.sign(disc.loc[observed,'mu'].to_numpy())
    agree=(np.sign(effects.loc[observed,'g'].to_numpy())==signs)
    pools={s:disc.index[(np.sign(disc.mu)==s)&(disc.index.isin(common))&(~disc.index.isin(selected.index))] for s in (-1,1)}
    need={s:int(sum(signs==s)) for s in (-1,1)}
    assert all(len(pools[s])>=need[s] for s in (-1,1))
    rng=np.random.default_rng(20260925); null=np.zeros(10000)
    for b in range(10000):
        picks=np.concatenate([rng.choice(pools[s],size=need[s],replace=False) for s in (-1,1)])
        null[b]=np.mean(np.sign(effects.loc[picks,'g'].to_numpy())==np.sign(disc.loc[picks,'mu'].to_numpy()))
    pd.DataFrame({'matched_random_fraction_agree':null}).to_csv(f'results/ic_urine_p7_null_{stratum}.csv.gz',index=False)
    score=float(agree.mean()); p=(1+int((null>=score).sum()))/(len(null)+1)
    rows.append(dict(stratum=stratum,n_cases=int(cases.sum()),n_controls=int(normal.sum()),n_selected=50,n_measured=len(observed),n_agree=int(agree.sum()),fraction_agree=score,null_mean=float(null.mean()),null_sd=float(null.std(ddof=1)),empirical_p=p,source=URL))
    for symbol,a in zip(observed,agree):
        gene_rows.append(dict(stratum=stratum,gene=symbol,discovery_mu=float(disc.loc[symbol,'mu']),new_g=float(effects.loc[symbol,'g']),new_v=float(effects.loc[symbol,'v']),agrees=int(a)))
os.makedirs('results',exist_ok=True)
pd.DataFrame(rows).to_csv('results/ic_urine_p7.csv',index=False)
pd.DataFrame(gene_rows).to_csv('results/ic_urine_p7_genes.csv',index=False)
pd.DataFrame({'gsm':ann.index,'title':ann.Sample_title.to_numpy(),'disease_state':state.to_numpy(),'group':np.where(normal,'control',np.where(without,'PBS_without_lesions','PBS_with_lesions'))}).to_csv('results/ic_urine_p7_labels.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
