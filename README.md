If you use this, cite GOrilla:

- Eran Eden*, Roy Navon*, Israel Steinfeld, Doron Lipson and Zohar Yakhini. "GOrilla: A Tool For Discovery And Visualization of Enriched GO Terms in Ranked Gene Lists", BMC Bioinformatics 2009, 10:48.
- Eran Eden, Doron Lipson, Sivan Yogev, Zohar Yakhini. "Discovering Motifs in Ranked Lists of DNA sequences", PLoS Computational Biology, 3(3):e39, 2007. 

## Usage

```python
from pygorilla import submit_ranked_genes

genes = ["TP53", "BRCA1", "EGFR", "MYC", "AKT1"]
result = submit_ranked_genes(genes, processing_timeout=30, poll_interval=1.0)

print(result.run_id)
for record in result.go_terms[:5]:
    print(record.go_id, record.description, record.p_value)
```

If the list is short or yields no enriched terms, the result will contain an empty `go_terms`
list and any warning messages that GOrilla returned under `result.metadata.warnings`.

### Target vs. background mode

Provide a background list to switch to the hypergeometric (two unranked list) mode automatically:

```python
target = ["TP53", "BRCA1", "EGFR", "MYC", "AKT1"]
background = target + ["PTEN", "KRAS", "CD4", "AKT2", "MAPK8"]

result = submit_ranked_genes(target, background_genes=background, analysis_name="demo run")
print(len(result.go_terms))  # -> number of enriched terms (possibly zero)
```

`fast_mode` applies only to the ranked-list workflow; it is ignored when a background list
is supplied.

### Working with pandas

Install pandas (already listed as a dependency) to export the enrichment results as DataFrames:

```python
normalized = result.to_go_terms_dataframe()
denormalized = result.to_go_terms_dataframe(normalized=False)
```

The normalized frame contains one row per GO term along with the aggregate counts and the lists
of gene symbols/descriptions. The denormalized version expands each GO term so that every gene
appears on its own row, making it convenient to join with other gene-level data.
