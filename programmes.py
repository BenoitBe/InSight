"""Programmes de réponse : les protéines significatives dans au moins un contraste
sont regroupées par la forme de leur profil sur les conditions (moyennes par
condition, centrées-réduites par protéine). Indépendant de toute notion de
temps : un programme est un groupe de protéines qui varient ensemble.

Chaque programme est ensuite décrit par :
  - les conditions où il est haut / bas ;
  - les thèmes fonctionnels sur-représentés parmi ses protéines (test
    hypergéométrique, univers = protéines significatives qui portent au moins
    un terme, BH q < 0.05 comme les GO de Deimos) ;
  - ses protéines les plus marquées portant un nom lisible.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

THEME_Q_THR = 0.05      # même seuil que les GO de Deimos (gProfiler / STRING, q < 5 %)


@dataclass
class Programmes:
    z: pd.DataFrame            # protéines x conditions (z-scores)
    labels: pd.Series          # protéine -> programme (1..k, ordonnés)
    k: int
    themes: dict               # programme -> [(nom du thème, p, n)]
    names: dict                # programme -> [noms lisibles]
    high: dict                 # programme -> condition(s) où il est le plus haut
    low: dict


def _choose_k(Z, kmin=3, kmax=8, seed=42):
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    best, best_s = kmin, -1
    n = len(Z)
    for k in range(kmin, min(kmax, max(kmin, n // 15)) + 1):
        lab = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(Z)
        s = silhouette_score(Z, lab, sample_size=min(n, 2000), random_state=seed)
        if s > best_s + 0.01:          # un k plus grand doit vraiment mieux séparer
            best, best_s = k, s
    return best


def build_programmes(run, design, themes=None, k=None, seed=42) -> Programmes | None:
    from sklearn.cluster import KMeans
    if run.means is None or len(run.conditions) < 3:
        return None
    sig = np.zeros(len(run.de), bool)
    for c in run.contrasts:
        sig |= run.sig(c).fillna(False).values
    ids = [p for p in run.de.index[sig] if p in run.means.index]
    if len(ids) < 20:
        return None
    conds = [c for c in design.conditions if c in run.means.columns]
    if len(conds) < 3:
        return None
    P = run.means.loc[ids, conds]
    sd = P.std(axis=1).replace(0, np.nan)
    Z = P.sub(P.mean(axis=1), axis=0).div(sd, axis=0).dropna()
    k = k or _choose_k(Z.values, seed=seed)
    km = KMeans(n_clusters=k, n_init=30, random_state=seed).fit(Z.values)
    raw = pd.Series(km.labels_, index=Z.index)
    # ordre : programmes les plus peuplés d'abord
    order = raw.value_counts().index.tolist()
    lab = raw.map({c: i + 1 for i, c in enumerate(order)})

    lfc_cols = [c + "_diff" for c in run.contrasts]
    peak = run.de.loc[Z.index, lfc_cols].abs().max(axis=1)
    names, high, low, th = {}, {}, {}, {}
    from .labels import informative
    for p in range(1, k + 1):
        members = lab.index[lab == p]
        prof = Z.loc[members].median()
        high[p] = prof.idxmax(); low[p] = prof.idxmin()
        cand = peak.loc[members].sort_values(ascending=False).index
        nm = [run.labels.get(i, i) for i in cand if informative(run.labels.get(i, ""))]
        names[p] = list(dict.fromkeys(nm))[:6]
    if themes is not None:
        th = _theme_enrichment(lab, themes)
    return Programmes(z=Z, labels=lab, k=k, themes=th, names=names, high=high, low=low)


def _theme_enrichment(lab, themes, max_themes=3):
    """Thèmes sur-représentés dans chaque programme (hypergéométrique + BH)."""
    from scipy.stats import hypergeom
    universe = set(lab.index)
    tg = {t: themes.genes(t) & universe for t in themes.names}
    tg = {t: g for t, g in tg.items() if len(g) >= 3}
    annotated = set().union(*tg.values()) if tg else set()
    if not annotated:
        return {}
    N = len(annotated)
    out = {}
    for p in sorted(lab.unique()):
        members = set(lab.index[lab == p]) & annotated
        n = len(members)
        rows = []
        for t, g in tg.items():
            k_ = len(g & members)
            if k_ >= 2:
                rows.append((t, hypergeom.sf(k_ - 1, N, len(g), n), k_))
        if rows:
            ps = np.array([r[1] for r in rows])
            order = np.argsort(ps)
            q = np.minimum.accumulate((ps[order] * len(ps) / (np.arange(len(ps)) + 1))[::-1])[::-1]
            sel = [(themes.names[rows[i][0]], float(qq), rows[i][2])
                   for i, qq in zip(order, q) if qq < THEME_Q_THR]
            out[p] = sel[:max_themes]
        else:
            out[p] = []
    return out
