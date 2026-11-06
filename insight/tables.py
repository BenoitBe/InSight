"""Tables d'InSight : ce que montrent les figures, en données réutilisables.

insight_tables.xlsx
  README              provenance (version, source, règle de sélection, seuils, k) et
                      dictionnaire des colonnes
  Facts               indicateurs et constats chiffrés du rapport, section par section
  Programmes          une ligne par protéine des programmes : programme, ajustement au
                      profil du programme (r, silhouette), thèmes, contrastes significatifs,
                      part de valeurs observées, moyennes log2 et z-scores par condition
  Programme_summary   une ligne par programme : taille, conditions haute / basse,
                      thèmes enrichis (q), protéines exemples
  Themes              une ligne par thème : termes, sources, protéines, contrastes et
                      programmes où il est enrichi
  Theme_terms         une ligne par terme : thème, source, FDR, protéines
  Theme_x_contrast    format long thème × contraste (FDR du meilleur terme, direction)

insight_programmes.tsv : la feuille Programmes en texte tabulé (lecture directe sous R).
Rien n'est recalculé hors des programmes et thèmes déjà construits pour les figures ;
r et silhouette décrivent la qualité d'affectation de chaque protéine.
"""
from __future__ import annotations

import html
import os
import re

import numpy as np
import pandas as pd

from . import __version__

XLSX_NAME = "insight_tables.xlsx"
TSV_NAME = "insight_programmes.tsv"
_TAG = re.compile(r"<[^>]+>")


def _plain(s: str) -> str:
    """Texte d'un constat du rapport, sans balises HTML."""
    return html.unescape(_TAG.sub("", str(s))).strip()


def _theme_id(t) -> str:
    return f"T{int(t) + 1}"


def _join(xs, sep="; ") -> str:
    return sep.join(str(x) for x in xs)


# --------------------------------------------------------------------------- programmes

def programmes_table(run, design, prog, themes=None) -> pd.DataFrame:
    """Une ligne par protéine affectée à un programme."""
    from sklearn.metrics import silhouette_samples
    Z = prog.z
    lab = prog.labels.reindex(Z.index)
    centroids = {p: Z[lab == p].mean() for p in range(1, prog.k + 1)}
    r = pd.Series({pid: float(np.corrcoef(Z.loc[pid].values, centroids[lab[pid]].values)[0, 1])
                   for pid in Z.index})
    sil = pd.Series(silhouette_samples(Z.values, lab.values), index=Z.index) \
        if prog.k > 1 and len(Z) > prog.k else pd.Series(np.nan, index=Z.index)

    # thèmes : ceux du programme (enrichis) qui contiennent la protéine, et tous les thèmes
    th_genes = {t: themes.genes(t) for t in themes.names} if themes is not None else {}
    enriched = {p: [name for name, _q, _n in prog.themes.get(p, [])] for p in range(1, prog.k + 1)}
    name_to_t = {v: t for t, v in themes.names.items()} if themes is not None else {}

    sig_cols = {c: run.sig(c).fillna(False) for c in run.contrasts}
    rows = []
    for pid in Z.index:
        p = int(lab[pid])
        mine = [t for t, g in th_genes.items() if pid in g]
        in_prog = [n for n in enriched[p] if pid in th_genes.get(name_to_t.get(n), set())]
        sig_in = [run.short(c) for c in run.contrasts if bool(sig_cols[c].get(pid, False))]
        lfc = {c: run.de.at[pid, c + "_diff"] for c in run.contrasts if pid in run.de.index}
        peak_c = max(lfc, key=lambda c: abs(lfc[c]) if np.isfinite(lfc[c]) else -1) if lfc else None
        rows.append(dict(
            protein=pid, name=run.labels.get(pid, pid) if run.labels is not None else pid,
            programme=f"P{p}", r_to_programme=round(r[pid], 3),
            silhouette=round(float(sil[pid]), 3) if np.isfinite(sil[pid]) else np.nan,
            high_in=prog.high[p], low_in=prog.low[p],
            programme_themes=_join(in_prog),
            all_themes=_join(themes.names[t] for t in mine) if mine else "",
            n_significant=len(sig_in), significant_in=_join(sig_in),
            peak_contrast=run.short(peak_c) if peak_c else "",
            peak_log2FC=round(float(lfc[peak_c]), 3) if peak_c else np.nan))
    out = pd.DataFrame(rows).set_index("protein")

    # part de valeurs observées (DPC) : un profil porté par des non-détections est signalé
    if run.observed is not None and run.samples is not None:
        obs = run.observed.reindex(Z.index)
        cond = run.samples.set_index("label")["condition"]
        cols = [c for c in obs.columns if c in cond.index]
        frac = obs[cols].T.groupby(cond.loc[cols].values).mean().T
        frac = frac.reindex(columns=[c for c in Z.columns if c in frac.columns])
        out["observed_fraction"] = obs[cols].mean(axis=1).round(3)
        out["min_observed_fraction"] = frac.min(axis=1).round(3)
        out["least_observed_in"] = frac.idxmin(axis=1)
    conds = list(Z.columns)
    means = run.means.reindex(index=Z.index, columns=conds)
    out = out.join(means.round(3).add_prefix("log2_"))
    out = out.join(Z.round(3).add_prefix("z_"))
    order = {f"P{p}": p for p in range(1, prog.k + 1)}
    out = (out.assign(_o=out.programme.map(order))
              .sort_values(["_o", "r_to_programme"], ascending=[True, False])
              .drop(columns="_o"))
    return out.reset_index()


def programme_summary(prog, ptab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for p in range(1, prog.k + 1):
        sub = ptab[ptab.programme == f"P{p}"]
        th = prog.themes.get(p, [])
        rows.append(dict(
            programme=f"P{p}", n_proteins=len(sub), high_in=prog.high[p], low_in=prog.low[p],
            mean_r_to_programme=round(float(sub.r_to_programme.mean()), 3),
            mean_silhouette=round(float(sub.silhouette.mean()), 3)
            if sub.silhouette.notna().any() else np.nan,
            enriched_themes=_join(f"{n} (q={q:.2g}, {k} prot.)" for n, q, k in th) or "none",
            examples=_join(prog.names.get(p, []), ", ")))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- thèmes

def themes_tables(run, themes, prog=None, labels=None):
    """(Themes, Theme_terms, Theme_x_contrast)."""
    terms = themes.terms.copy()
    size = terms.groupby("theme").term.nunique()
    best = themes.table.groupby("theme").p.min()
    in_prog = {}
    if prog is not None:
        for p, lst in prog.themes.items():
            for n, q, _k in lst:
                in_prog.setdefault(n, []).append(f"P{p} (q={q:.2g})")
    lab = labels if labels is not None else {}
    rows = []
    for t in sorted(themes.names, key=lambda t: best.get(t, 1.0)):
        sub = terms[terms.theme == t].sort_values("p")
        genes = sorted(themes.genes(t))
        tc = themes.table[themes.table.theme == t].sort_values("p")
        rows.append(dict(
            theme_id=_theme_id(t), theme=themes.names[t], n_terms=int(size.get(t, 0)),
            n_proteins=len(genes), best_fdr=float(best.get(t, np.nan)),
            sources=_join(sorted(set(sub.source)), ", "),
            enriched_in=_join(f"{run.short(c)} (FDR={p:.2g})" for c, p in zip(tc.contrast, tc.p)),
            programmes=_join(in_prog.get(themes.names[t], [])),
            terms=_join(f"{i} {n}" for i, n in zip(sub.term, sub.name)),
            proteins=_join((lab.get(g, g) for g in genes), ", "),
            protein_ids=_join(genes, ", ")))
    tt = terms.assign(theme_id=terms.theme.map(_theme_id),
                      theme_name=terms.theme.map(themes.names),
                      n_proteins=terms.genes.map(len),
                      proteins=terms.genes.map(lambda g: _join(sorted(g), ", ")))
    tt = tt.sort_values(["theme", "p"])[["theme_id", "theme_name", "term", "name", "source",
                                          "p", "n_proteins", "proteins"]]
    tt = tt.rename(columns={"term": "term_id", "name": "term_name", "p": "best_fdr"})
    tx = themes.table.assign(theme_id=themes.table.theme.map(_theme_id),
                             theme_name=themes.table.theme.map(themes.names),
                             comparison=themes.table.contrast.map(run.short))
    tx = tx.sort_values(["theme", "p"])[["theme_id", "theme_name", "contrast", "comparison",
                                          "p", "z", "n"]]
    tx = tx.rename(columns={"p": "best_fdr", "z": "direction_z", "n": "n_terms"})
    return pd.DataFrame(rows), tt.reset_index(drop=True), tx.reset_index(drop=True)


# --------------------------------------------------------------------------- écriture

_DICT = [
    ("Programmes", "protein", "identifiant Deimos (tête du Protein.Group)"),
    ("Programmes", "name", "nom lisible (symbole, orthologue, description ou accession)"),
    ("Programmes", "programme", "programme de réponse (k-means sur les z-scores des moyennes)"),
    ("Programmes", "r_to_programme", "corrélation de Pearson entre le profil z de la protéine et "
                                     "le profil moyen de son programme"),
    ("Programmes", "silhouette", "silhouette de la protéine (-1 à 1) : < 0 = plus proche d'un "
                                 "autre programme"),
    ("Programmes", "high_in / low_in", "conditions où le programme est le plus haut / bas"),
    ("Programmes", "programme_themes", "thèmes enrichis dans le programme qui contiennent la "
                                       "protéine"),
    ("Programmes", "all_themes", "tous les thèmes dont un terme contient la protéine"),
    ("Programmes", "n_significant / significant_in", "comparaisons où la protéine passe la "
                                                     "règle de sélection de Deimos"),
    ("Programmes", "peak_contrast / peak_log2FC", "comparaison au |log2FC| maximal"),
    ("Programmes", "observed_fraction", "part des échantillons où la protéine est détectée (DPC)"),
    ("Programmes", "min_observed_fraction / least_observed_in",
     "part observée minimale sur les conditions, et laquelle : un profil porté par des "
     "non-détections (valeurs estimées) y apparaît"),
    ("Programmes", "log2_<condition>", "moyenne log2 par condition (matrice de l'abondance)"),
    ("Programmes", "z_<condition>", "z-score de la moyenne (profil utilisé par le k-means)"),
    ("Themes", "theme_id", "T1 = thème au meilleur terme ; mêmes thèmes que la figure 04"),
    ("Themes", "best_fdr", "FDR du meilleur terme du thème (tables GO de Deimos, q < 0.05)"),
    ("Themes", "enriched_in", "comparaisons où au moins un terme du thème est enrichi"),
    ("Themes", "programmes", "programmes où le thème est sur-représenté (q < 0.05)"),
    ("Theme_x_contrast", "direction_z", "z-score médian des termes (> 0 : protéines plutôt "
                                        "en hausse)"),
]


def _format(writer):
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    head = PatternFill("solid", fgColor="F0EFEC")
    for ws in writer.book.worksheets:
        if ws.max_row < 1:
            continue
        heads = [1] + ([r for r in range(2, ws.max_row + 1) if ws.cell(r, 1).value == "sheet"]
                       if ws.title == "README" else [])
        for r in heads:                                  # en-têtes (dont le dictionnaire)
            for cell in ws[r]:
                cell.font = Font(bold=True)
                cell.fill = head
                cell.alignment = Alignment(vertical="top", wrap_text=False)
        ws.freeze_panes = "B2" if ws.title in ("Programmes", "Themes", "Theme_terms") else "A2"
        if ws.max_row > 1 and ws.title not in ("README",):
            ws.auto_filter.ref = ws.dimensions
        for j, col in enumerate(ws.iter_cols(min_row=1, max_row=min(ws.max_row, 200)), 1):
            width = max((len(str(c.value)) for c in col if c.value is not None), default=8)
            ws.column_dimensions[get_column_letter(j)].width = min(max(9, width + 2), 60)


def write_tables(out_dir, run, design, prog=None, themes=None, sections=(), kpis=(),
                 meta=(), k_requested=None) -> str | None:
    """Écrit insight_tables.xlsx (+ insight_programmes.tsv). Renvoie le chemin du classeur."""
    sheets = {}
    labels = run.labels.to_dict() if run.labels is not None else {}

    facts = [dict(section="Key figures", fact=f"{v} {n}") for n, v in kpis]
    for s in sections:
        facts += [dict(section=s["title"], fact=_plain(b)) for b in s.get("bullets", [])]
    sheets["Facts"] = pd.DataFrame(facts, columns=["section", "fact"])

    ptab = None
    if prog is not None:
        ptab = programmes_table(run, design, prog, themes)
        sheets["Programmes"] = ptab
        sheets["Programme_summary"] = programme_summary(prog, ptab)
    if themes is not None:
        th, tt, tx = themes_tables(run, themes, prog, labels)
        sheets["Themes"], sheets["Theme_terms"], sheets["Theme_x_contrast"] = th, tt, tx

    from .io import GO_Q_THR
    from .programmes import THEME_Q_THR
    from .themes import JACCARD_MIN
    prov = [("InSight version", __version__),
            ("Source workbook", os.path.basename(run.xlsx)),
            ("Deimos output folder", os.path.abspath(run.out_dir)),
            ("Protein selection", f"{run.rule_text()} [{run.rule_source}]"),
            ("GO / pathway terms", f"q < {GO_Q_THR:g} (tables GO de Deimos)"),
            ("Themes", f"terms grouped when Jaccard >= {JACCARD_MIN:g} (shared proteins), "
                       f"modularity communities"),
            ("Abundance", run.abundance_kind or "n/a"),
            ("Protein names", run.label_source)]
    if prog is not None:
        prov += [("Programmes", f"k = {prog.k} ({'set by user' if k_requested else 'silhouette, 3-8'}"
                                f"); k-means on z-scored condition means of the "
                                f"{len(prog.z)} selected proteins"),
                 ("Programme themes", f"hypergeometric in the selected proteins, BH q < "
                                      f"{THEME_Q_THR:g}")]
    prov += [(m, "") for m in meta]
    readme = pd.DataFrame(prov, columns=["item", "value"])
    dico = pd.DataFrame(_DICT, columns=["sheet", "column", "meaning"])

    path = os.path.join(out_dir, XLSX_NAME)
    with pd.ExcelWriter(path, engine="openpyxl") as w:
        readme.to_excel(w, sheet_name="README", index=False)
        dico.to_excel(w, sheet_name="README", index=False, startrow=len(readme) + 2)
        for name, df in sheets.items():
            df.to_excel(w, sheet_name=name, index=False)
        _format(w)
    if ptab is not None:
        ptab.to_csv(os.path.join(out_dir, TSV_NAME), sep="\t", index=False)
    return path
