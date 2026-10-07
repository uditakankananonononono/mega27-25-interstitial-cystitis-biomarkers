# Source provenance and rights documentation gaps (interstitial-cystitis child, 2026-10-08 audit record)

This file records what this repository's own documents say about source hashes and terms. It is a documentation record, not a license verdict, not an integrity certificate, and not a clearance of any recorded rights hold. A matching claim below means two recorded hash strings are equal at a named field in pinned repository documents. It does not mean any assay bytes were re-downloaded or verified. No rights record was found for these sources, so reuse and redistribution are NOT cleared by this file.

Audited child commit: `54697488b2162488ec30ac1a300896dc46170de0`  
Compared shared-parent repo: `mega27-25-biomarkers-underserved-diseases` at `36ff87995f4fe0e1d08f9e2b2cd885d341b88937`

Rights: No source-specific rights check found in inspected disease-child documentation. Hash matches are not legal clearance or byte verification.

## 1. Result-file hash records in this child (2)

Hash values are shown as 12-hex display prefixes only; the full value is in the cited file and field. Classes are kept separate on purpose: an expression/assay source hash, a GEO series-matrix hash (metadata that may include expression), a reference annotation, and a published artifact are different kinds of record.

- `results/ic_p25_result.json` `/matrix_sha256` hash `cdc2ca391b68...` (GSE238208); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/ic_p25_result.json` `/expression_sha256` hash `c4b8c08bbc41...` (GSE238208); class: expression_or_assay_source_hash_record; equal to shared-parent field

## 4. Metadata crosswalk tables (not assay bytes)

These per-sample tables carry hash columns describing GEO source-response metadata. They are not blanket expression-matrix integrity.

- `projects/interstitial_cystitis/sources/GSE57560_used_sample_crosswalk.csv`: columns ['sha256'], 7 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/interstitial_cystitis/sources/GSE238208_used_sample_crosswalk.csv`: columns ['sha256'], 50 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/interstitial_cystitis/sources/GSE11783_used_sample_crosswalk.csv`: columns ['sha256'], 16 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/interstitial_cystitis/sources/GSE621_historically_used_sample_crosswalk.csv`: columns ['sha256'], 14 records; file identical to a shared-parent file: True; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/interstitial_cystitis/sources/gse621_treatment/GSE621_treatment_crosswalk.csv`: columns ['sha256', 'source_matrix_sha256'], 8 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity
- `projects/interstitial_cystitis/sources/gse57560_normal_capacity/GSE57560_normal_capacity_used_crosswalk.csv`: columns ['source_sha256', 'matrix_sha256'], 9 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity

## Coverage boundary

- Source accessions listed in the child manifest: 5.
- Of those, 4 have no mapped result-file payload hash in this audit: GSE11783, GSE28242, GSE57560, GSE621.
- 0 hash fields inherited from other diseases' records are not counted toward this child.
- Records absent from child manifest can still have result-file hashes. Neither presence nor absence proves assay acquisition by this scout. Other-disease copied hashes never count toward this child.
