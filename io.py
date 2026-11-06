"""Lecture d'un run Deimos terminé (dossier de sortie) : résultats différentiels,
contrastes, règle volcano, matrices, scores PCA, tables GO, orthologues.

Tout est optionnel sauf la feuille Differential_Expression : chaque figure
vérifie ce dont elle a besoin et se désactive proprement sinon.
"""
from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class Run:
    out_dir: str
    xlsx: str
    de: pd.DataFrame                       # une ligne par protéine, index = protein id
    contrasts: list                        # noms de contrastes, ordre du manifest
    pairs: dict                            # contraste -> (numérateur, dénominateur)
    p_thr: float = 0.05
    ratio_min: float = 1.5                 # défauts de Deimos ; la règle du run prime
    use_padj: bool = False
    rule_source: str = "Deimos defaults"
    y_axis: str = "pvalue"                 # volcano_y_axis de Deimos
    samples: pd.DataFrame | None = None    # label, condition, (facteurs)
    conditions: list = field(default_factory=list)   # ordre d'affichage
    means: pd.DataFrame | None = None      # protéines x conditions (log2)
    abundance: pd.DataFrame | None = None  # protéines x échantillons (log2)
    observed: pd.DataFrame | None = None   # protéines x échantillons (1/0)
    abundance_kind: str = ""
    scores: pd.DataFrame | None = None     # PCA_Scores
    go: dict = field(default_factory=dict)  # contraste -> table GO
    labels: pd.Series | None = None        # protein id -> nom lisible
    label_source: str = ""
    config: dict = field(default_factory=dict)

    @property
    def lfc_min(self) -> float:
        return float(np.log2(self.ratio_min)) if self.ratio_min > 0 else 0.0

    def sig(self, c: str) -> pd.Series:
        """Protéines significatives pour un contraste : la règle volcano de Deimos
        (volcano_status), statistique < seuil ET |log2FC| >= log2(ratio)."""
        p = self.de[c + ("_p.adj" if self.use_padj else "_p.val")]
        return (p < self.p_thr) & (self.de[c + "_diff"].abs() >= self.lfc_min)

    def rule_text(self) -> str:
        return (f"{'p.adj' if self.use_padj else 'p.value'} < {self.p_thr:g} & "
                f"|log2FC| >= {self.lfc_min:.3g} (ratio >= {self.ratio_min:g})")

    def short(self, c: str) -> str:
        num, den = self.pairs.get(c, (None, None))
        return f"{num} vs {den}" if num else c.replace("_vs_", " vs ")


def _find_xlsx(out_dir: str) -> str:
    cands = sorted(glob.glob(os.path.join(out_dir, "ProteomicAnalysis_Results*.xlsx")),
                   key=os.path.getmtime)
    if not cands:
        raise FileNotFoundError(f"No ProteomicAnalysis_Results*.xlsx in {out_dir}")
    return cands[-1]


def _manifest(xl: pd.ExcelFile) -> pd.DataFrame | None:
    if "_manifest" not in xl.sheet_names:
        return None
    m = pd.read_excel(xl, sheet_name="_manifest", header=None).iloc[1:, :3]
    m.columns = ["section", "key", "value"]
    return m


# Termes GO / voies : q < 0.05, le seuil de l'enrichissement de Deimos (gProfiler, STRING)
GO_Q_THR = 0.05


def _volcano_rule(text: str):
    """'p.adj < 0.05 & ratio >= 2.0' -> (use_padj, p, ratio)."""
    t = str(text)
    use_padj = "p.adj" in t or "fdr" in t.lower()
    mp = re.search(r"<\s*([0-9.eE-]+)", t)
    mr = re.search(r"ratio\s*>=\s*([0-9.]+)", t)
    return use_padj, float(mp.group(1)) if mp else 0.05, float(mr.group(1)) if mr else 1.0


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _rule(man, cfg):
    """Règle de sélection de Deimos : clés structurées du manifeste (schéma 1.2), sinon
    son texte volcano_rule, sinon last_config.yaml, sinon les défauts de Deimos."""
    st = {} if man is None else {r.key: r.value for r in man.itertuples()
                                 if r.section == "stats"}
    if _num(st.get("volcano_p_thresh")) is not None:
        return (str(st.get("volcano_p_column", "p.adj")).strip() == "p.adj",
                _num(st["volcano_p_thresh"]), _num(st.get("volcano_ratio_min")) or 1.0,
                "Deimos manifest")
    if st.get("volcano_rule"):
        return (*_volcano_rule(st["volcano_rule"]), "Deimos manifest")
    if cfg and _num(cfg.get("volcano_p_thresh")) is not None:
        return (bool(cfg.get("volcano_use_padj", False)), _num(cfg["volcano_p_thresh"]),
                _num(cfg.get("volcano_ratio_min")) or 1.5, "last_config.yaml")
    return False, 0.05, 1.5, "Deimos defaults (no rule found)"


def _read_tsv(path: str) -> pd.DataFrame | None:
    return pd.read_csv(path, sep="\t") if os.path.exists(path) else None


def _config(out_dir: str) -> dict:
    p = os.path.join(out_dir, "last_config.yaml")
    if not os.path.exists(p):
        return {}
    try:
        import yaml
        with open(p, encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}
    except Exception:                                                 # noqa: BLE001
        return {}


def load_run(out_dir: str, xlsx: str | None = None, reference_fasta: str | None = None,
             verbose: bool = True) -> Run:
    xlsx = xlsx or _find_xlsx(out_dir)
    xl = pd.ExcelFile(xlsx)
    de = pd.read_excel(xl, sheet_name="Differential_Expression").copy()
    id_col = "name" if "name" in de.columns else de.columns[0]
    de["protein"] = de[id_col].astype(str)
    de = de.drop_duplicates("protein").set_index("protein", drop=False)

    man = _manifest(xl)
    pairs, contrasts, go_map = {}, [], {}
    cfg = _config(out_dir)
    use_padj, p_thr, ratio, rule_src = _rule(man, cfg)
    if man is not None:
        for r in man.itertuples():
            if r.section == "contrast" and isinstance(r.value, str) and "|" in r.value:
                pairs[r.key] = tuple(r.value.split("|", 1))
                contrasts.append(r.key)
            elif r.section == "go_sheet":
                go_map[r.key] = r.value
    if not contrasts:                     # sans manifest : colonnes *_diff / *_p.adj
        contrasts = [c[:-5] for c in de.columns if c.endswith("_diff")
                     and not c.endswith("_imputed_diff") and c[:-5] + "_p.adj" in de.columns]
        pairs = {c: tuple(c.split("_vs_", 1)) for c in contrasts if "_vs_" in c}
    contrasts = [c for c in contrasts if c + "_diff" in de.columns]

    mdir = os.path.join(out_dir, "matrices")
    samples = _read_tsv(os.path.join(mdir, "samples.tsv"))
    run = Run(out_dir=out_dir, xlsx=xlsx, de=de, contrasts=contrasts, pairs=pairs,
              p_thr=p_thr, ratio_min=ratio, use_padj=use_padj, config=cfg,
              rule_source=rule_src,
              y_axis=str(cfg.get("volcano_y_axis", "pvalue")).lower() if cfg else "pvalue")

    # --- échantillons / conditions ---
    if samples is None and "PCA_Scores" in xl.sheet_names:
        s = pd.read_excel(xl, sheet_name="PCA_Scores")
        samples = s[["label", "condition"]].copy()
    if samples is not None:
        samples["label"] = samples["label"].astype(str)
        samples["condition"] = samples["condition"].astype(str)
        run.samples = samples
        run.conditions = list(dict.fromkeys(samples["condition"]))

    # --- matrices : DPC si disponible, sinon log2 imputée ---
    ab = _read_tsv(os.path.join(mdir, "dpc_abundance.tsv"))
    kind = "DPC"
    if ab is None:
        ab, kind = _read_tsv(os.path.join(mdir, "log2_imputed.tsv")), "log2 imputed"
    if ab is None and "Log2_Impute" in xl.sheet_names and samples is not None:
        li = pd.read_excel(xl, sheet_name="Log2_Impute")
        li["protein"] = li["Protein.Group"].astype(str).str.split(";").str[0]
        ab, kind = li, "log2 imputed"
    if ab is not None and samples is not None:
        pid = "protein" if "protein" in ab.columns else ab.columns[0]
        cols = [c for c in samples["label"] if c in ab.columns]
        A = ab.set_index(ab[pid].astype(str))[cols].astype(float)
        A = A[~A.index.duplicated()]
        run.abundance, run.abundance_kind = A, kind
        ob = _read_tsv(os.path.join(mdir, "dpc_observed.tsv")) if kind == "DPC" else None
        if ob is not None:
            run.observed = ob.set_index(ob["protein"].astype(str))[cols].astype(int)
        cond = samples.set_index("label").loc[cols, "condition"]
        run.means = A.T.groupby(cond.values).mean().T.reindex(columns=run.conditions)
        run.means = run.means.loc[:, run.means.notna().any()]

    if "PCA_Scores" in xl.sheet_names:
        run.scores = pd.read_excel(xl, sheet_name="PCA_Scores")

    # --- tables GO, rattachées aux contrastes ---
    for sheet, c in go_map.items():
        if sheet in xl.sheet_names and c in run.contrasts:
            t = pd.read_excel(xl, sheet_name=sheet)
            t = t.loc[:, ~t.columns.astype(str).str.startswith("Unnamed")]
            if "term_id" in t.columns and len(t.dropna(subset=["term_id"])):
                t = t.dropna(subset=["term_id"])
                if "p_value" in t.columns:              # termes significatifs (q < 5 %)
                    t = t[pd.to_numeric(t["p_value"], errors="coerce") < GO_Q_THR]
                if len(t):
                    run.go[c] = t.reset_index(drop=True)

    from .labels import best_labels
    run.labels, run.label_source = best_labels(run, xl, reference_fasta)
    if verbose:
        print(f"  [InSight] {len(de)} proteins, {len(run.contrasts)} contrasts, "
              f"{len(run.conditions)} conditions, GO tables: {len(run.go)}, "
              f"abundance: {run.abundance_kind or 'none'}, names: {run.label_source}")
        print(f"  [InSight] protein selection: {run.rule_text()} [{run.rule_source}]; "
              f"GO terms: q < {GO_Q_THR:g}")
    return run
