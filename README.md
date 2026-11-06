# InSight — biological reading of a Deimos run

InSight turns a finished [Deimos](https://github.com/BenoitBe/deimos-proteomics-pipeline) run into
a short visual story: what changes, where the conditions sit, which groups of proteins move
together, which functions they carry, and the named proteins behind each pattern. It reads the
results workbook and the `matrices/` folder and **recomputes no statistics**: the differential
calls, enrichments and abundances are exactly those of the Deimos run.

It works with **any design**. Time courses get trajectories and line profiles; designs with no
order (treatments, tissues, genotypes, sexes) get labelled condition maps and bar profiles. An order
is read from condition names only when there is evidence for it (see below), never assumed.

## Output

Figures (PNG + SVG for editing) and a self-contained `insight_report.html` (images embedded, no
network, light/dark), with factual bullet points computed from the results next to each figure.

| # | Figure | What it shows |
|---|--------|---------------|
| 01 | `01_landscape` | Up / down counts per comparison, with the run's volcano rule |
| 02 | `02_condition_map` | PCA centroids (replicates faint, 1-SD ellipses). Ordered design with ≤ 6 groups: arrows in level order. More groups: centroids coloured by level. No order: labelled centroids |
| 03 | `03_programmes` | *Response programmes*: differential proteins clustered by the shape of their profile across conditions (k-means on z-scored condition means, k by silhouette). Each programme is described by enriched themes (hypergeometric test within the differential proteins, BH q < 0.05). Layout: lines per group (ordered, ≤ 6 groups), bars per condition (no order, ≤ 14 conditions), heatmap otherwise |
| 04 | `04_themes` | *Functional themes*: GO/KEGG/Reactome terms from all comparisons merged when they share proteins (Jaccard ≥ 0.35, modularity communities). A: enrichment map; B: theme × comparison matrix (size = best FDR, colour = direction z-score) |
| 05 | `05_volcano_NN_*` | Volcanos of the comparisons with most differential proteins (or the ones you choose), coloured by the top themes, with the most marked named proteins |
| 06 | `06_spotlight` | One or two named proteins per programme, replicate by replicate (hollow dots = not detected, DPC estimate) |

**Thresholds are those of Deimos.** The proteins are the ones Deimos selected: the volcano rule
of the run (statistic `p.adj` or raw `p.value`, threshold, ratio: `p < thr & |log2FC| >=
log2(ratio)`), read from the `_manifest` (structured keys of schema 1.2, else its `volcano_rule`
text), else from `last_config.yaml`, else Deimos' defaults (raw p < 0.05, ratio ≥ 1.5). GO /
pathway terms are kept at **q < 0.05**, the gProfiler / STRING threshold; the themes enriched in a
programme use the same 5 %. Volcano Y axis and threshold line follow Deimos' `volcano_y_axis`.

Protein names: informative gene symbol → orthologue symbol (Perseverance + `GN=` of the reference
FASTA, or STRING) → short description → accession. Locus tags, `LOC…`, UniProt accessions and
RefSeq IDs are never used as names when something better exists.

## Usage

```bash
python insight.py proteomics_output                       # Deimos output folder
python insight.py proteomics_output --title "P5 — B. napus imbibition"
python insight.py proteomics_output --volcano A_vs_B --volcano C_vs_D --k 6
python insight.py proteomics_output --categorical         # never draw trajectories
python insight.py proteomics_output --order-regex "T(\d+)"
python insight.py proteomics_output --reference-fasta reference.fasta
```

From Deimos (last step, after Pathfinder, separate process, never blocks the run):
`make_insight: true` in the YAML, or `python deimos.py --insight`.

| Option | Default | |
|---|---|---|
| `--out` | `<deimos_out>/insight` | output folder |
| `--xlsx` | latest `ProteomicAnalysis_Results*.xlsx` | workbook to read |
| `--reference-fasta` | auto (Perseverance keys of `last_config.yaml`, then the run and data folders) | reference proteome for orthologue names |
| `--k` | silhouette over 3–8 | number of response programmes |
| `--volcano` | the `--max-volcanos` (4) comparisons with most differential proteins | repeatable |
| `--categorical` | off | ignore any order in condition names |
| `--order-regex` | automatic | regex whose first group is the ordered level |

## How the design is read

A number in a condition name becomes an ordered level only with evidence:

- a unit stuck to it: `24h`, `7d`, `48hai`, `15DAP`, `30min`, `2wk` …
- an explicit prefix: `T3`, `D7`, `J14`, `Day7`, `Dose10`, `Stage2` …
- or a crossed design: at least two groups sharing at least two levels (`embAviso0/2/8/24` and
  `embMajor0/2/8/24`).

There must also be ≥ 3 distinct levels and ≥ 2 levels per group. Otherwise the conditions stay
categorical: `Line1/Line2/Line3`, `Grp_Male4/Grp_Female4`, `S1…S4`, `WT_clone2` are identifiers,
not an order. Common prefixes/suffixes are dropped for display (`Tissue_Liver` → `Liver`).

## What it needs

Only the `Differential_Expression` sheet is mandatory; every figure checks its inputs and is skipped
with a message otherwise. Richer runs give more: `_manifest` (contrasts, volcano rule, GO sheet
mapping), `GO_*` sheets (any backend: gProfiler, STRING, local Perseverance annotation),
`PCA_Scores`, `matrices/` (`samples.tsv`, `dpc_abundance.tsv` + `dpc_observed.tsv`, else
`log2_imputed.tsv`), `Orthologues_Perseverance` / `Orthologues_STRING`.

Python ≥ 3.10 with `numpy pandas scipy scikit-learn matplotlib networkx adjustText openpyxl
pyyaml` (`pip install -r requirements.txt`).

## Reading the figures — limits

- Programmes are clusters of profile *shapes* among differential proteins; k is a choice (forced
  with `--k`). Two programmes can be close; read them with the themes and the spotlight.
- Themes merge redundant terms; a theme named after its most significant term may contain related
  terms (see the enrichment map). Direction comes from the z-score computed by Deimos.
- Bullet points are facts computed from the results. The biological interpretation is yours.

## Tests

```bash
python -m pytest -q tests
```

Synthetic Deimos outputs (`tests/fixture.py`) cover categorical designs, 18 conditions, time
courses, crossed designs with many groups, two conditions and a workbook reduced to one sheet.
