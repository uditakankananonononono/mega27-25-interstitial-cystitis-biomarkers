# Record-unit derivation (2026-10-08)

Documentation note, not a license verdict and not an integrity certificate. It explains the "disease-tagged record units" count in README.md.

## What a record unit is
One row of `results/disease_tagged_accession_manifest.csv` (columns: accession, source, role). A row is one GEO accession record that this project used or fetched for a stated role: an individual GSM sample record or a GSE series record. Rows are not patients, not independent datasets and not independent cohorts. Many GSM rows are nested in a parent GSE and the `role` text says so. No per-patient or per-cell unit is counted.

## Current count
123 rows, 123 unique accessions: {'other': 1, 'GSE series': 5, 'GSM sample': 117}. Non-GSE/GSM row: EFO_1000869.

## Why README said 106
The README count was written when the manifest had 106 rows and was not updated when later commits added rows:
manifest rows by commit: 9665123=106, acdc124=115, 7c3544a=123.
- acdc124 (2026-09-28) +9 GSM rows: exploratory low- vs normal-capacity IC probe-effect comparison, nested in GSE57560.
- 7c3544a (2026-09-28) +8 GSM rows: exploratory APF vs mock in a cultured normal bladder cell line (GSE621 source-verified assays).
Those commits changed the manifest, not the README sentence, so the README lagged by 17. The earlier 106 equals the manifest at commit 9665123.

## Method
`git show <commit>:results/disease_tagged_accession_manifest.csv`, row count and (accession, role) set difference between commits; no row was removed between commits.

## Boundary
The recount does not change what the rows mean. The README caveat stays: these are nested record units, not independent datasets or people.
