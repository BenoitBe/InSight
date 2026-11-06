"""Thèmes fonctionnels : les termes d'enrichissement (GO, KEGG, Reactome...) de
tous les contrastes sont regroupés quand ils reposent sur les mêmes protéines
(indice de Jaccard), par détection de communautés. Un thème = un groupe de
termes redondants, nommé par son terme le plus significatif (GO:BP de
préférence)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

JACCARD_MIN = 0.35
MIN_GENES = 3


@dataclass
class Themes:
    terms: pd.DataFrame          # term, name, source, p, genes (set), theme
    hits: pd.DataFrame           # contrast, term, p, z, genes, theme
    names: dict                  # theme -> nom
    graph: object = None         # networkx.Graph des termes
    table: pd.DataFrame = field(default_factory=pd.DataFrame)  # theme x contrast

    def top(self, n=16, min_terms=1) -> list:
        best = self.table.groupby("theme").p.min().sort_values()
        size = self.terms.groupby("theme").term.nunique()
        return [t for t in best.index if size.get(t, 0) >= min_terms][:n]

    def genes(self, theme) -> set:
        sub = self.terms[self.terms.theme == theme]
        return set().union(*sub.genes) if len(sub) else set()


def _genes(x) -> set:
    if isinstance(x, (list, tuple, set)):
        return {str(g).strip() for g in x if str(g).strip()}
    return {g.strip() for g in str(x).replace(";", ",").split(",") if g.strip()
            and g.strip().lower() != "nan"}


def build_themes(go: dict, jaccard_min: float = JACCARD_MIN) -> Themes | None:
    import networkx as nx
    rows = []
    for c, t in go.items():
        icol = "intersection" if "intersection" in t.columns else (
            "intersections" if "intersections" in t.columns else None)
        if icol is None:
            continue
        for r in t.itertuples(index=False):
            d = r._asdict()
            rows.append(dict(contrast=c, term=str(d["term_id"]),
                             name=str(d.get("term_name") or d["term_id"]),
                             source=str(d.get("source", "")),
                             p=float(d.get("p_value", np.nan)),
                             z=float(d.get("z_score", 0) or 0),
                             genes=_genes(d[icol])))
    if not rows:
        return None
    hits = pd.DataFrame(rows)
    terms = (hits.groupby("term")
             .agg(name=("name", "first"), source=("source", "first"), p=("p", "min"),
                  genes=("genes", lambda s: set().union(*s)))
             .reset_index())
    terms = terms[terms.genes.map(len) >= MIN_GENES].reset_index(drop=True)
    if terms.empty:
        return None
    G = nx.Graph()
    G.add_nodes_from(terms.term)
    gs = dict(zip(terms.term, terms.genes))
    keys = list(gs)
    for i, a in enumerate(keys):
        ga = gs[a]
        for b in keys[i + 1:]:
            gb = gs[b]
            inter = len(ga & gb)
            if inter:
                j = inter / len(ga | gb)
                if j >= jaccard_min:
                    G.add_edge(a, b, w=j)
    comms = nx.algorithms.community.greedy_modularity_communities(G, weight="w")
    theme_of, names = {}, {}
    info = terms.set_index("term")
    for k, cset in enumerate(sorted(comms, key=lambda c: info.loc[list(c)].p.min())):
        sub = info.loc[list(cset)]
        pick = sub[sub.source == "GO:BP"] if (sub.source == "GO:BP").any() else sub
        names[k] = pick.sort_values(["p", "name"]).name.iloc[0]
        for t in cset:
            theme_of[t] = k
    terms["theme"] = terms.term.map(theme_of)
    hits["theme"] = hits.term.map(theme_of)
    hits = hits.dropna(subset=["theme"])
    hits["theme"] = hits.theme.astype(int)
    table = (hits.groupby(["theme", "contrast"])
             .agg(p=("p", "min"), z=("z", "median"), n=("term", "nunique"))
             .reset_index())
    return Themes(terms=terms, hits=hits, names=names, graph=G, table=table)
