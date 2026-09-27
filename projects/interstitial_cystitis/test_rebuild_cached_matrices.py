from rebuild_cached_matrices import rebuild

def test_GSE621_rebuild():
    r=rebuild('GSE621')
    assert r['raw_probe_by_GSM']==[4132,28]
    assert r['cached_gene_by_GSM']==[3556,14]
    assert r['max_absolute_value_delta']==0

def test_GSE11783_rebuild():
    r=rebuild('GSE11783')
    assert r['raw_probe_by_GSM']==[54675,16]
    assert r['cached_gene_by_GSM']==[21594,16]
    assert r['max_absolute_value_delta']<1e-6
