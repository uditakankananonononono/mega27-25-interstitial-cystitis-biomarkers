import csv
import gzip
import hashlib
from pathlib import Path
from replay import run

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def test_source_and_treatment_unit():
    r=run()
    assert (r['n_used_GSM'],r['n_cell_lines'],r['n_APF_libraries'],r['n_mock_libraries'])==(8,2,4,4)
    assert r['n_probes']==r['finite_probe_pairs']==4132
    cross=list(csv.DictReader(open(HERE/'GSE621_treatment_crosswalk.csv')))
    assert len({v['gsm'] for v in cross})==8
    assert {v['cell_line'] for v in cross}=={'1','2'}
    for row in cross:
        raw=gzip.open(HERE/f'{row["gsm"]}.soft.txt.gz','rb').read()
        assert hashlib.sha256(raw).hexdigest()==row['sha256']
    manifest=list(csv.DictReader(open(ROOT/'results/disease_tagged_accession_manifest.csv')))
    assert len(manifest)==123 and len({v['accession'] for v in manifest})==123
    assert all('not clinical IC diagnosis' in next(x['role'] for x in manifest if x['accession']==v['gsm']) for v in cross)
