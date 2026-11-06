# Changelog

## 0.1.1 — 2026-11-06

Thresholds aligned on Deimos.
- Protein selection read from the Deimos manifest: structured keys (schema 1.2) first, then the
  `volcano_rule` text, then `last_config.yaml`, then Deimos' defaults (raw p < 0.05, ratio ≥ 1.5;
  previously p.adj < 0.05, ratio ≥ 2 when no rule was found). Rule and source printed and
  written in the report.
- GO / pathway terms kept at q < 0.05 (rows above are dropped, whatever the backend).
- Themes enriched in a programme: BH q < 0.05 (was 0.1).
- Volcano Y axis and threshold line placed as in Deimos (`volcano_y_axis`, raw-p equivalent of
  the p.adj threshold).

## 0.1.0 — 2026-11-06

First version: biological reading of a finished Deimos run, any design.

- Reads the Deimos workbook (`Differential_Expression`, `_manifest`, `GO_*`, `PCA_Scores`,
  orthologue sheets) and `matrices/`; nothing is recomputed. Only `Differential_Expression` is
  mandatory, every figure is skipped cleanly when its inputs are missing.
- Design detection from condition names with evidence only (unit, explicit prefix, or crossed
  design); categorical otherwise. Overrides: `--categorical`, `--order-regex`, `order_column`.
- Figures 01–06 (landscape, condition map, response programmes, functional themes, themed
  volcanos, spotlight), PNG + SVG, and a self-contained HTML report with computed facts.
- Readable protein names: gene symbol, then orthologue symbol (Perseverance + reference FASTA, or
  STRING), then description, then accession.
- Layouts adapt to the design: trajectories / lines for ordered designs with ≤ 6 groups, level
  colours beyond; bars for categorical designs, heatmap above 14 conditions.
- Tests on synthetic outputs (categorical, 18 conditions, time course, crossed design with 8
  groups, 2 conditions, single-sheet workbook) and on design detection.
