"""P25 registered lesion/non-lesion paired IC signature stress test."""
import csv,hashlib,json,re,sys,urllib.request
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests
from ubiomark.geo import parse_series_matrix
root=Path(__file__).resolve().parents[1]
folder=root/'data/geo/p25';folder.mkdir(parents=True,exist_ok=True)
gse='GSE238208';stem='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238208/'
sources=[(folder/(gse+'_series_matrix.txt.gz'),stem+'matrix/GSE238208_series_matrix.txt.gz','cdc2ca391b688abebe47d5ab4e8f89e3d0800612bac04ae704f3a4fda2ed45e5'),(folder/'GSE238208_gene_expression_yoshiyukiakiyama.txt.gz',stem+'suppl/GSE238208_gene_expression_yoshiyukiakiyama.txt.gz','c4b8c08bbc418b30e800b289e3ab28ffaf6a49b071440ce092e2dfbe9d49cc4a')]
for p,url,sha in sources:
 if not p.exists():urllib.request.urlretrieve(url,p)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
_,ann,_=parse_series_matrix(str(sources[0][0]));assert len(ann)==63 and ann.index.is_unique
prior={r['accession'] for r in csv.DictReader((root/'results/dataset_manifest.csv').open())};assert not prior.intersection(ann.index)
name_to_gsm=dict(zip(ann.Sample_description,ann.index));assert len(name_to_gsm)==63
x=pd.read_csv(sources[1][0],sep='\t',skiprows=1,low_memory=False)
fpkm=[c for c in x if c.endswith('.FPKM')]
assert len(fpkm)==63 and {c[:-5] for c in fpkm}==set(name_to_gsm)
vals=x[fpkm].to_numpy(dtype=float);assert np.isfinite(vals).all() and (vals>=0).all()
assert len(x)>10000
symbols=x['Name'].fillna('').astype(str)
valid=(symbols.str.fullmatch(r'[A-Za-z][A-Za-z0-9._-]*') & ~symbols.duplicated(keep=False))
expr=pd.DataFrame(np.log2(vals[valid]+1),index=symbols[valid],columns=[name_to_gsm[c[:-5]] for c in fpkm]);assert expr.index.is_unique
pair=[]
for i in range(1,26):
 a=ann.index[ann.Sample_title==f'Hunner_Type_Interstitial_Cystitis_Case{i}_Hunner_Lesion'].tolist()
 b=ann.index[ann.Sample_title==f'Hunner_Type_Interstitial_Cystitis_Case{i}_Non_Hunner_Lesion'].tolist()
 assert len(a)==len(b)==1
 assert ann.loc[a[0],'Sample_characteristics_ch1_2']=='biopsied site: Hunner_Lesion'
 assert ann.loc[b[0],'Sample_characteristics_ch1_2']=='biopsied site: Non_Hunner_Lesion'
 assert all('Hunner_Type_Interstitial_Cystitis' in ann.loc[z,'Sample_characteristics_ch1_1'] for z in (a[0],b[0]))
 pair.append(dict(patient=i,lesion_gsm=a[0],nonlesion_gsm=b[0],lesion_column=ann.loc[a[0],'Sample_description']+'.FPKM',nonlesion_column=ann.loc[b[0],'Sample_description']+'.FPKM'))
assert len({r[k] for r in pair for k in ('lesion_gsm','nonlesion_gsm')})==50
ids1=[r['lesion_gsm'] for r in pair];ids0=[r['nonlesion_gsm'] for r in pair]
other=set(ann.index)-set(ids1)-set(ids0);assert len(other)==13 and all('BCG' in ann.loc[z,'Sample_title'] for z in other)
disc=pd.read_csv(root/'results/meta_discovery/interstitial_cystitis.csv.gz',index_col=0)
selected=disc.sort_values('p').head(50)
difference=expr[ids1].to_numpy()-expr[ids0].to_numpy()
effects=pd.Series(difference.mean(axis=1),index=expr.index)
common=disc.index.intersection(expr.index)
observed=selected.index.intersection(common)
need={sign:int((np.sign(disc.loc[observed,'mu'])==sign).sum()) for sign in [-1,1]}
assert all(disc.loc[observed,'k']>=selected.k.min())
match=np.sign(effects.loc[observed])==np.sign(disc.loc[observed,'mu'])
pools={sign:disc.index[(disc.index.isin(common))&(~disc.index.isin(selected.index))&(np.sign(disc.mu)==sign)&(disc.k>=selected.k.min())].to_numpy() for sign in [-1,1]}
assert all(len(pools[s])>=need[s] for s in [-1,1])
rng=np.random.default_rng(20260925)
null=np.zeros(10000,dtype=float)
for i in range(len(null)):
 picks=np.concatenate([rng.choice(pools[s],size=need[s],replace=False) for s in [-1,1]])
 null[i]=np.mean(np.sign(effects.loc[picks])==np.sign(disc.loc[picks,'mu']))
score=float(match.mean());p=float((1+(null>=score).sum())/(len(null)+1))
t=pd.Series({g:float(ttest_rel(expr.loc[g,ids1],expr.loc[g,ids0]).pvalue) for g in observed})
q=multipletests(t.values,method='fdr_bh')[1]
rows=[dict(gene=g,discovery_mu=float(disc.loc[g,'mu']),discovery_k=int(disc.loc[g,'k']),paired_mean_log2_fpkm1_diff=float(effects[g]),agrees=int(match.loc[g]),paired_t_p=float(t[g]),panel_BH_q=float(q[i])) for i,g in enumerate(observed)]
with (root/'results/ic_p25_sample_map.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(pair[0]));w.writeheader();w.writerows(pair)
with (root/'results/ic_p25_genes.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
pd.DataFrame({'random_fraction_agree':null}).to_csv(root/'results/ic_p25_null.csv.gz',index=False)
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238208',matrix_sha256=sources[0][2],expression_sha256=sources[1][2],n_paired_people=len(pair),n_bcg_excluded=len(other),n_unique_gene_symbols=len(expr),n_panel_measured=len(observed),n_sign_match=int(match.sum()),fraction_agree=score,null_mean=float(null.mean()),empirical_p=p,registered_descriptive_rule_pass=bool(len(observed)>=40 and int(match.sum())>=40 and p<.01),n_panel_BH_q_under_05=int((q<.05).sum()),panel_nominal_p_under_05=int((t<.05).sum()))
(root/'results/ic_p25_result.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
