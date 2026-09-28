import csv
import gzip
import hashlib
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path.insert(0,str(HERE))
from fetch_and_analyze import run,PIN


def test_capacity_source_and_results():
    result=run()
    assert (result['n_normal_capacity_IC'],result['n_low_capacity_IC'],result['n_healthy_controls_excluded'])==(9,4,3)
    assert result['n_probes']==62976 and result['n_finite_effects']==62976
    rows=list(csv.DictReader(open(HERE/'GSE57560_normal_capacity_used_crosswalk.csv')))
    assert len({r['gsm'] for r in rows})==9
    for r in rows:
        raw=gzip.open(HERE/f'{r["gsm"]}.soft.txt.gz','rb').read()
        assert hashlib.sha256(raw).hexdigest()==r['source_sha256']
        assert r['matrix_sha256']==PIN
    manifests=[list(csv.DictReader(open(ROOT/'results'/n))) for n in ['dataset_manifest.csv','disease_tagged_accession_manifest.csv']]
    assert len(manifests[0])==len(manifests[1])==115
    assert {r['accession'] for r in manifests[0]}=={r['accession'] for r in manifests[1]}
    assert {r['gsm'] for r in rows}.issubset({r['accession'] for r in manifests[0]})
