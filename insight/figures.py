"""Figures InSight. Chaque fonction renvoie (chemins écrits, faits pour le texte)
ou None si les données nécessaires manquent. Style sobre et commun ; palette
catégorielle à ordre fixe (identité), divergente bleu-gris-rouge (sens)."""
from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID, FAINT = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#d6d5d0"
DIV = LinearSegmentedColormap.from_list("insight_div", ["#2a78d6", "#f0efec", "#e34948"])
SEQ = LinearSegmentedColormap.from_list("insight_seq", ["#e3eefb", "#2a78d6", "#0b2f5c"])
UP, DOWN = "#e34948", "#2a78d6"


def style():
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11,
        "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.labelcolor": INK2, "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
        "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False,
        "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
        "figure.facecolor": "white", "axes.facecolor": "white", "svg.fonttype": "none",
    })


def _save(fig, out_dir, name, formats=("png", "svg")):
    paths = []
    for ext in formats:
        p = os.path.join(out_dir, f"{name}.{ext}")
        fig.savefig(p)
        paths.append(p)
    plt.close(fig)
    return paths


def _title(ax, title, subtitle=None, extra=0):
    """Titre + sous-titre au-dessus d'un axe ; extra = points réservés entre l'axe
    et le sous-titre (étiquettes posées au-dessus du cadre)."""
    from matplotlib.transforms import offset_copy
    ax.set_title(title, pad=(18 if subtitle else 6) + extra)
    if subtitle:
        tr = offset_copy(ax.transAxes, fig=ax.figure, y=extra, units="points")
        ax.text(0, 1.01, subtitle, transform=tr, fontsize=8, color=INK2, va="bottom")


def _suptitle(fig, title, subtitle=None, y=1.0):
    """Titre + sous-titre au-dessus des panneaux ; écarts fixés en pouces pour ne
    jamais se chevaucher, quelle que soit la hauteur de la figure."""
    h = fig.get_figheight()
    y_sub = y + 0.03 / h
    y_tit = y_sub + (0.20 / h if subtitle else 0)
    fig.text(0.01, y_tit, title, fontsize=13, fontweight="bold", color=INK, va="bottom")
    if subtitle:
        fig.text(0.01, y_sub, subtitle, fontsize=8, color=INK2, va="bottom")


def _short(s, n=34):
    s = str(s)
    return s if len(s) <= n else s[:n - 1] + "…"


def group_colors(design):
    """Couleur d'identité par groupe, ordre fixe ; au-delà de 8 groupes, pas de
    couleur générée : une seule teinte, les groupes sont nommés en clair."""
    if not design.group:
        return {}
    gs = design.groups
    if len(gs) > len(CAT):
        return {g: CAT[0] for g in gs}
    return {g: CAT[i] for i, g in enumerate(gs)}


# ---------------------------------------------------------------------------
# 1. Paysage des contrastes
# ---------------------------------------------------------------------------
def fig_landscape(run, out_dir):
    rows = []
    for c in run.contrasts:
        s = run.sig(c).fillna(False)
        lf = run.de[c + "_diff"]
        rows.append((c, int((s & (lf > 0)).sum()), int((s & (lf < 0)).sum())))
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=["contrast", "up", "down"])
    h = max(2.4, 0.32 * len(df) + 1.2)
    fig, ax = plt.subplots(figsize=(7.6, h))
    y = np.arange(len(df))[::-1]
    ax.barh(y, df.up, color=UP, height=0.62, label="up")
    ax.barh(y, -df.down, color=DOWN, height=0.62, label="down")
    m = max(df.up.max(), df.down.max(), 1)
    for yi, u, d in zip(y, df.up, df.down):
        ax.text(u + m * 0.015, yi, str(u), va="center", fontsize=7.5, color=INK2)
        ax.text(-d - m * 0.015, yi, str(d), va="center", ha="right", fontsize=7.5, color=INK2)
    ax.set_yticks(y, [_short(run.short(c), 40) for c in df.contrast], fontsize=8)
    ax.set_xlim(-m * 1.18, m * 1.18)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xticks(ax.get_xticks(), [f"{abs(int(t))}" for t in ax.get_xticks()])
    ax.set_xlabel("significant proteins (down | up)")
    ax.grid(axis="y", visible=False)
    rule = (f"{'FDR' if run.use_padj else 'p'} < {run.p_thr:g}, ratio ≥ {run.ratio_min:g}")
    _title(ax, "How much changes in each comparison", f"Differential proteins ({rule})")
    ax.legend(loc="lower right", fontsize=8)
    top = df.assign(tot=df.up + df.down).sort_values("tot", ascending=False)
    facts = {"n_contrasts": len(df), "empty": int(((df.up + df.down) == 0).sum()),
             "top": [(run.short(r.contrast), r.up, r.down) for r in top.head(3).itertuples()]}
    return _save(fig, out_dir, "01_landscape"), facts


# ---------------------------------------------------------------------------
# 2. Carte des conditions (PCA) — trajectoires si le plan est ordonné
# ---------------------------------------------------------------------------
def fig_map(run, design, out_dir):
    from matplotlib.patches import Ellipse, FancyArrowPatch
    sc = run.scores
    if sc is None or "PC1" not in sc.columns:
        return None
    sc = sc.copy()
    sc["condition"] = sc["condition"].astype(str)
    v1 = sc.filter(like="PC1_var").iloc[0, 0] if len(sc.filter(like="PC1_var").columns) else None
    v2 = sc.filter(like="PC2_var").iloc[0, 0] if len(sc.filter(like="PC2_var").columns) else None
    cen = sc.groupby("condition")[["PC1", "PC2"]].mean()
    gcol = group_colors(design)
    many = len(design.groups) > 3
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    span = np.ptp(sc[["PC1", "PC2"]].values, axis=0).max() or 1
    for c in design.conditions:
        if c not in cen.index:
            continue
        col = gcol.get(design.group[c], CAT[0]) if design.group else CAT[0]
        r = sc[sc.condition == c]
        ax.scatter(r.PC1, r.PC2, s=22, color=col, alpha=0.35, lw=0, zorder=2)
        if len(r) >= 3:
            w, V = np.linalg.eigh(np.cov(r[["PC1", "PC2"]].T))
            ang = np.degrees(np.arctan2(V[1, 1], V[0, 1]))
            ax.add_patch(Ellipse(cen.loc[c].values, 2 * np.sqrt(max(w[1], 0)),
                                 2 * np.sqrt(max(w[0], 0)), angle=ang, fc=col, ec="none",
                                 alpha=0.08, zorder=1))
    facts = {"var": (v1, v2)}
    traj = design.ordered and len(design.groups or ["all"]) <= 6
    if traj:
        groups = design.groups or ["all"]
        steps = []
        for g in groups:
            cs = [c for c in design.conditions if (design.group or {}).get(c, "all") == g
                  and c in cen.index]
            cs = sorted(cs, key=lambda c: design.level[c])
            col = gcol.get(g, CAT[0])
            for a, b in zip(cs[:-1], cs[1:]):
                ax.add_patch(FancyArrowPatch(cen.loc[a].values, cen.loc[b].values,
                                             arrowstyle="-|>", mutation_scale=12, lw=1.8,
                                             color=col, shrinkA=6, shrinkB=6, zorder=3))
                steps.append((g, a, b, float(np.linalg.norm(cen.loc[b] - cen.loc[a]))))
            sizes = np.linspace(35, 90, len(cs))
            ax.scatter(cen.loc[cs].PC1, cen.loc[cs].PC2, s=sizes, color=col, ec="white",
                       lw=1.4, zorder=4)
            # étiquettes de niveau ; niveaux confondus regroupés
            grp = [[cs[0]]] if cs else []
            for c in cs[1:]:
                if np.linalg.norm(cen.loc[c] - cen.loc[grp[-1][-1]]) < 0.08 * span:
                    grp[-1].append(c)
                else:
                    grp.append([c])
            for gg in grp:
                x, y = cen.loc[gg].mean()
                ax.annotate(", ".join(design.label(design.level[c]) for c in gg), (x, y),
                            xytext=(0, 11), textcoords="offset points", ha="center",
                            fontsize=8, color=INK2, zorder=5,
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                                      alpha=0.8))
            if cs and design.group:
                # nom du groupe dans le prolongement du dernier pas (jamais sur la flèche)
                x, y = cen.loc[cs[-1]]
                if len(cs) > 1:
                    d = (cen.loc[cs[-1]] - cen.loc[cs[-2]]).values / [np.ptp(sc.PC1) or 1,
                                                                   np.ptp(sc.PC2) or 1]
                    d = d / (np.linalg.norm(d) or 1)
                else:
                    d = np.array([1.0, 0.0])
                ax.annotate(_short(g, 22), (x, y), xytext=(14 * d[0], 14 * d[1] - 4),
                            textcoords="offset points",
                            ha="left" if d[0] >= -0.3 else "right",
                            va="bottom" if d[1] >= 0 else "top",
                            fontsize=10, fontweight="bold", color=col, zorder=5)
        facts["steps"] = sorted(steps, key=lambda s: -s[3])[:3]
        sub = (f"Condition centroids linked in order of {design.level_name} · "
               "replicates faint · ellipses = 1 SD")
    else:
        # sans trajectoire : centroïdes nommés ; si le plan est ordonné mais a trop de
        # groupes pour des flèches lisibles, la couleur porte le niveau (séquentiel)
        texts = []
        present = [c for c in design.conditions if c in cen.index]
        if design.ordered:
            lv = design.levels
            cols = {l: SEQ(0.25 + 0.75 * i / max(1, len(lv) - 1)) for i, l in enumerate(lv)}
            for l in lv:
                ax.scatter([], [], s=50, color=cols[l], ec="white", label=design.label(l))
            ax.legend(title=design.level_name, fontsize=7.5, title_fontsize=8, loc="best")
        if len(present) <= 16:
            named = set(present)
        else:   # trop de noms : les 12 conditions les plus excentrées
            dist = np.linalg.norm(cen.loc[present] - cen.loc[present].mean(), axis=1)
            named = set(np.array(present)[np.argsort(-dist)[:12]])
        for c in present:
            col = cols[design.level[c]] if design.ordered else CAT[0]
            ax.scatter(*cen.loc[c], s=70, color=col, ec="white", lw=1.4, zorder=4)
            if c in named:
                x, y = cen.loc[c]
                texts.append(ax.text(x + 0.012 * span, y + 0.012 * span,
                                     _short(design.display(c), 24), fontsize=8, color=INK,
                                     zorder=5))
        ax.margins(0.10)
        try:
            from adjustText import adjust_text
            cc = cen.reindex([c for c in design.conditions if c in cen.index])
            adjust_text(texts, x=cc.PC1.values, y=cc.PC2.values, ax=ax, expand=(1.3, 1.7),
                        force_static=(0.3, 0.5), iter_lim=500,
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.4))
        except Exception:                                             # noqa: BLE001
            pass
        sub = ("Condition centroids · replicates faint · ellipses = 1 SD"
               + (f" · colour = {design.level_name} ({len(design.groups)} groups, too many "
                  "for trajectories)" if design.ordered else "")
               + (" · 12 most distinct named" if len(present) > 16 else ""))
    d = cen.values
    if len(d) >= 2:
        D = np.linalg.norm(d[:, None] - d[None], axis=2)
        i, j = np.unravel_index(D.argmax(), D.shape)
        facts["farthest"] = (cen.index[i], cen.index[j])
    ax.set_xlabel(f"PC1 ({v1:.1f} %)" if v1 is not None else "PC1")
    ax.set_ylabel(f"PC2 ({v2:.1f} %)" if v2 is not None else "PC2")
    ax.axhline(0, color=GRID, lw=0.8, zorder=0); ax.axvline(0, color=GRID, lw=0.8, zorder=0)
    ax.margins(0.14)
    _title(ax, "Where the conditions sit in protein space", sub)
    return _save(fig, out_dir, "02_condition_map"), facts


# ---------------------------------------------------------------------------
# 3. Programmes de réponse
# ---------------------------------------------------------------------------
def fig_programmes(run, design, prog, out_dir):
    if prog is None:
        return None
    k = prog.k
    conds = [c for c in design.conditions if c in prog.z.columns]
    gcol = group_colors(design)
    lines = design.ordered and len(design.levels) >= 3 and len(design.groups or ["all"]) <= 6
    if not lines and len(conds) > 14:
        return _programme_heatmap(run, design, prog, out_dir)
    ncol = 3 if k > 4 else 2
    nrow = int(np.ceil(k / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.0 * nrow + 0.4), sharey=True,
                             squeeze=False)
    for p in range(1, k + 1):
        ax = axes.flat[p - 1]
        Z = prog.z[prog.labels == p]
        if lines:
            lv = design.levels
            x = np.arange(len(lv))
            for g in (design.groups or ["all"]):
                cs = {design.level[c]: c for c in conds
                      if (design.group or {}).get(c, "all") == g}
                cols = [cs[l] for l in lv if l in cs]
                xs = [x[lv.index(l)] for l in lv if l in cs]
                m = Z[cols]
                col = gcol.get(g, CAT[0])
                ax.fill_between(xs, m.quantile(0.25), m.quantile(0.75), color=col, alpha=0.15, lw=0)
                ax.plot(xs, m.median(), color=col, lw=2, marker="o", ms=4, mec="white",
                        mew=0.8, label=_short(g, 18) if design.group else None)
            ax.set_xticks(x, [design.label(l) for l in lv])
        else:
            # conditions sans ordre : barres depuis 0 (pas de ligne qui suggérerait un ordre)
            x = np.arange(len(conds))
            med, q1, q3 = Z[conds].median(), Z[conds].quantile(0.25), Z[conds].quantile(0.75)
            ax.bar(x, med, width=0.62, color=[UP if v > 0 else DOWN for v in med], alpha=0.85,
                   zorder=2)
            ax.vlines(x, q1, q3, color=INK2, lw=0.9, zorder=3)
            ax.set_xticks(x, [_short(design.display(c), 14) for c in conds], rotation=45,
                          ha="right", fontsize=7.5)
            ax.grid(axis="x", visible=False)
        ax.axhline(0, color=MUTED, lw=0.6, ls=":")
        th = prog.themes.get(p, [])
        line2 = ("; ".join(_short(t[0], 30) for t in th[:2]) if th
                 else ", ".join(prog.names.get(p, [])[:4]))
        ax.set_title(f"P{p} · {len(Z)} proteins", fontsize=10, pad=15)
        ax.text(0, 1.02, _short(line2, 58), transform=ax.transAxes, fontsize=7.3,
                color=INK2, va="bottom", style="normal" if th else "italic")
    for ax in axes.flat[k:]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("abundance (z-score)")
    if lines and design.group:
        axes.flat[0].legend(fontsize=7.5, loc="best", handlelength=1.3)
    _suptitle(fig, "Response programmes: proteins that move together",
              f"{len(prog.z)} differential proteins grouped by profile shape (k-means, k={k}) · "
              + ("median ± IQR" if lines else "bar = median, whisker = IQR")
              + " · subtitle = enriched themes, or proteins")
    fig.tight_layout(rect=[0, 0, 1, 0.985], h_pad=2.0)
    return _save(fig, out_dir, "03_programmes"), {"layout": "lines" if lines else "bars"}


def _programme_heatmap(run, design, prog, out_dir):
    conds = [c for c in design.conditions if c in prog.z.columns]
    M = pd.DataFrame({p: prog.z[prog.labels == p][conds].median() for p in range(1, prog.k + 1)}).T
    fig, ax = plt.subplots(figsize=(min(16, 2.5 + 0.28 * len(conds)), 0.55 * prog.k + 2.2))
    lim = max(1.0, np.nanmax(np.abs(M.values)))
    im = ax.imshow(M.values, aspect="auto", cmap=DIV, norm=TwoSlopeNorm(0, -lim, lim))
    ax.set_xticks(range(len(conds)), [_short(design.display(c), 16) for c in conds], rotation=60,
                  ha="right", fontsize=7)
    ylab = []
    for p in M.index:
        th = prog.themes.get(p, [])
        extra = _short(th[0][0], 34) if th else ", ".join(prog.names.get(p, [])[:3])
        ylab.append(f"P{p} ({int((prog.labels == p).sum())}) · {extra}")
    ax.set_yticks(range(len(M)), ylab, fontsize=8)
    ax.grid(False)
    if design.group:
        gs = [design.group[c] for c in conds]
        starts = [i for i in range(len(gs)) if i == 0 or gs[i] != gs[i - 1]]
        for i in starts[1:]:
            ax.axvline(i - 0.5, color="white", lw=2.5)
        for a, b in zip(starts, starts[1:] + [len(gs)]):
            ax.text((a + b - 1) / 2, -0.62, _short(gs[a], max(6, int(4.5 * (b - a)))), ha="center",
                    va="bottom", fontsize=7, color=INK2)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("median z-score", fontsize=8)
    _title(ax, "Response programmes: proteins that move together",
           f"{len(prog.z)} differential proteins grouped by profile shape (k={prog.k}) · "
           "cell = median z-score of the programme", extra=12 if design.group else 0)
    return _save(fig, out_dir, "03_programmes"), {"layout": "heatmap"}


# ---------------------------------------------------------------------------
# 4. Thèmes fonctionnels
# ---------------------------------------------------------------------------
def _packed_layout(H, seed=7, aspect=1.2, pad=1.1):
    """Disposition « carte d'enrichissement » : chaque composante connexe est
    disposée séparément (force-directed) puis les composantes sont rangées en
    étagères, les plus grandes d'abord. Évite les amas isolés éparpillés au hasard
    et les étiquettes qui tombent sur l'amas voisin."""
    import networkx as nx
    comps = sorted(nx.connected_components(H), key=len, reverse=True)
    blocks = []
    for cset in comps:
        sub = H.subgraph(cset)
        n = len(cset)
        r = 0.35 + 0.30 * np.sqrt(n)
        if n == 1:
            p = {next(iter(cset)): np.zeros(2)}
        else:
            p = nx.spring_layout(sub, k=1.2 / np.sqrt(n), iterations=300, seed=seed, weight="w")
            P = np.array(list(p.values()))
            P = P - P.mean(0)
            P = P / (np.abs(P).max() or 1) * r
            p = dict(zip(p.keys(), P))
        blocks.append((p, r))
    total = sum((2 * r + pad) ** 2 for _, r in blocks)
    width = max(2 * blocks[0][1] + pad, np.sqrt(total / aspect) * 1.1) if blocks else 1
    pos, x, y, row_h = {}, 0.0, 0.0, 0.0
    for p, r in blocks:
        d = 2 * r + pad
        if x > 0 and x + d > width:
            x, y, row_h = 0.0, y - row_h, 0.0
        for k_, v in p.items():
            pos[k_] = v + np.array([x + d / 2, y - d / 2])
        x += d
        row_h = max(row_h, d)
    return pos


def fig_themes(run, themes, out_dir, n_themes=16):
    import networkx as nx
    if themes is None:
        return None
    keep = themes.top(n_themes)
    if not keep:
        return None
    order = [c for c in run.contrasts if c in set(themes.table.contrast)]
    tab = themes.table[themes.table.theme.isin(keep)]
    mat = tab.pivot(index="theme", columns="contrast", values="p").reindex(columns=order)
    mat = -np.log10(mat).fillna(0)
    if len(mat) > 2:
        from scipy.cluster.hierarchy import leaves_list, linkage
        ordr = list(mat.index[leaves_list(linkage(mat.values, "average"))])
    else:
        ordr = list(mat.index)
    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
    fig = plt.figure(figsize=(13, max(5.6, 0.36 * len(ordr) + 2.4)))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.75 + 0.06 * len(order)], wspace=0.12)
    size_t = themes.terms.groupby("theme").term.nunique()

    # A. carte d'enrichissement
    axA = fig.add_subplot(gs[0])
    nodes = [t for t in themes.graph.nodes
             if t in set(themes.terms.term[themes.terms.theme.isin(keep)])]
    H = themes.graph.subgraph(nodes)
    pos = _packed_layout(H, aspect=max(0.6, (0.36 * len(ordr) + 1.4) / 5.6))
    gsz = dict(zip(themes.terms.term, themes.terms.genes.map(len)))
    th_of = dict(zip(themes.terms.term, themes.terms.theme))
    zt = themes.hits.groupby("term").z.median()
    for a, b, d in H.edges(data=True):
        axA.plot(*zip(pos[a], pos[b]), color="#c9c8c3", lw=0.4 + 1.6 * d["w"], zorder=1)
    xy = np.array([pos[n] for n in H.nodes])
    axA.scatter(xy[:, 0], xy[:, 1], s=[14 + 5 * gsz[n] for n in H.nodes],
                c=[zt.get(n, 0) for n in H.nodes], cmap=DIV, norm=norm, ec="white", lw=0.8,
                zorder=3)
    texts = []
    for t in keep:
        pts = np.array([pos[n] for n in H.nodes if th_of[n] == t])
        if len(pts):
            # médoïde : un nœud du thème (la moyenne peut tomber entre deux amas)
            med = pts[np.argmin(((pts[:, None] - pts[None]) ** 2).sum(-1).sum(1))]
            ylow = pts[:, 1].min() - 0.18
            texts.append(axA.text(med[0], ylow, _short(themes.names[t], 26), fontsize=7.3,
                                  ha="center", va="top", zorder=4,
                                  bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                                            alpha=0.85)))
    x0, x1 = xy[:, 0].min(), xy[:, 0].max(); y0, y1 = xy[:, 1].min(), xy[:, 1].max()
    axA.set_aspect("equal", adjustable="box", anchor="N")
    axA.set_xlim(x0 - 0.22 * (x1 - x0 or 1) - 0.5, x1 + 0.22 * (x1 - x0 or 1) + 0.5)
    axA.set_ylim(y0 - 0.08 * (y1 - y0 or 1) - 0.4, y1 + 0.08 * (y1 - y0 or 1) + 0.4)
    try:
        from adjustText import adjust_text
        adjust_text(texts, x=xy[:, 0], y=xy[:, 1], ax=axA, expand=(1.15, 1.4),
                    min_arrow_len=2,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    except Exception:                                                 # noqa: BLE001
        pass
    axA.set_axis_off()
    _title(axA, "A  Enrichment map",
           "terms linked by shared proteins · size = proteins · colour = direction")

    # B. matrice thème x contraste
    axB = fig.add_subplot(gs[1])
    for i, t in enumerate(ordr):
        for j, c in enumerate(order):
            r = tab[(tab.theme == t) & (tab.contrast == c)]
            if r.empty:
                continue
            lp = -np.log10(max(float(r.p.iloc[0]), 1e-300))
            axB.scatter(j, i, s=30 + 70 * min(lp, 6), c=[float(r.z.iloc[0])], cmap=DIV,
                        norm=norm, ec=INK2, lw=0.4, zorder=3)
    axB.set_yticks(range(len(ordr)),
                   [f"{_short(themes.names[t], 42)}  ({size_t.get(t, 0)})" for t in ordr],
                   fontsize=8)
    axB.yaxis.tick_right()
    axB.set_xticks(range(len(order)), [_short(run.short(c), 24) for c in order], rotation=45,
                   ha="right", fontsize=7.5)
    axB.set_xlim(-0.6, len(order) - 0.4); axB.set_ylim(-0.7, len(ordr) - 0.3)
    axB.spines[["left", "bottom"]].set_visible(False)
    _title(axB, "B  Functional themes across comparisons",
           "size = −log10(FDR) of the best term · colour = direction (z-score) · (n terms)")
    hs = [axB.scatter([], [], s=30 + 70 * v, c="white", ec=INK2, lw=0.6) for v in (2, 4, 6)]
    pa = axA.get_position()
    fig.legend(hs, ["1e-2", "1e-4", "1e-6"], title="best-term FDR (panel B)", ncol=3,
               loc="center left", bbox_to_anchor=(pa.x0 + 0.52 * pa.width, pa.y0 - 0.035),
               fontsize=7.5, title_fontsize=8, handletextpad=0.9, columnspacing=2.0,
               borderpad=0.8)
    cax = fig.add_axes([pa.x0, pa.y0 - 0.045, 0.40 * pa.width, 0.018])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=DIV, norm=norm), cax=cax, orientation="horizontal")
    cb.set_label("direction (z-score): down ← → up", fontsize=8)
    facts = {"n_themes": len(themes.names), "n_terms": len(themes.terms),
             "top": [(themes.names[t], [run.short(c) for c in
                                        tab[tab.theme == t].sort_values("p").contrast.head(3)])
                     for t in keep[:5]]}
    return _save(fig, out_dir, "04_themes"), facts


# ---------------------------------------------------------------------------
# 5. Volcans annotés
# ---------------------------------------------------------------------------
def fig_volcano(run, c, themes, out_dir, idx):
    from .labels import informative
    de = run.de
    lfc = de[c + "_diff"]
    # axe Y comme les volcanos de Deimos (volcano_y_axis) : p brute par défaut
    y_padj = run.use_padj and run.y_axis == "padj"
    p = de[c + ("_p.adj" if y_padj else "_p.val")].clip(lower=1e-300)
    lp = -np.log10(p)
    s = run.sig(c).fillna(False)
    up, dn = s & (lfc > 0), s & (lfc < 0)
    groups = []
    if themes is not None:
        t = themes.table[themes.table.contrast == c].sort_values("p")
        used = set()
        for th in t.theme:
            genes = themes.genes(th) & set(de.index[s])
            genes -= used
            if len(genes) >= 3:
                groups.append((_short(themes.names[th], 34), CAT[len(groups)],
                               de.index.isin(list(genes))))
                used |= genes
            if len(groups) == 3:
                break
    fig, ax = plt.subplots(figsize=(7.8, 6.0))
    ns = ~s
    ax.scatter(lfc[ns], lp[ns], s=5, color=FAINT, lw=0, zorder=1, rasterized=True,
               label=f"not significant ({int(ns.sum())})")
    if groups:
        other = s.copy()
        for _, _, m in groups:
            other &= ~m
        ax.scatter(lfc[other], lp[other], s=11, color=MUTED, lw=0, zorder=2,
                   label=f"other significant ({int(other.sum())})")
        for name, col, m in groups:
            ax.scatter(lfc[m], lp[m], s=28, color=col, ec="white", lw=0.6, zorder=3,
                       label=f"{name} ({int(m.sum())})")
    else:
        ax.scatter(lfc[up], lp[up], s=14, color=UP, lw=0, zorder=2, label=f"up ({int(up.sum())})")
        ax.scatter(lfc[dn], lp[dn], s=14, color=DOWN, lw=0, zorder=2,
                   label=f"down ({int(dn.sum())})")
    stat = de[c + ("_p.adj" if run.use_padj else "_p.val")]
    # ligne de seuil placée comme dans volcano_axis de Deimos
    cut = None
    if y_padj or not run.use_padj:
        cut = run.p_thr
    elif (stat < run.p_thr).any():
        hi = float(p[stat < run.p_thr].max())
        above = p[(stat >= run.p_thr) & (p > hi)]
        cut = float(np.sqrt(hi * above.min())) if len(above) else hi
    if cut is not None:
        ax.axhline(-np.log10(cut), color=MUTED, lw=0.8, ls="--")
        ax.text(0.995, -np.log10(cut), f"{'p.adj' if run.use_padj else 'p'} < {run.p_thr:g}",
                fontsize=7.5, color=INK2, ha="right", va="bottom",
                transform=ax.get_yaxis_transform())
    for v in (run.lfc_min, -run.lfc_min):
        if v:
            ax.axvline(v, color=MUTED, lw=0.8, ls="--")
    xm = float(np.nanmax(np.abs(lfc.values))) if len(lfc) else 1.0
    ax.set_xlim(-xm * 1.22, xm * 1.22)
    ax.set_ylim(-0.02 * float(np.nanmax(lp)), float(np.nanmax(lp)) * 1.08)
    named = pd.Series([informative(run.labels.get(i, "")) for i in de.index], index=de.index)
    score = lp * lfc.abs()
    cand = list(score[s & named].sort_values(ascending=False).index[:10])
    for _, _, m in groups:
        cand += [i for i in score[m & named].sort_values(ascending=False).index[:3]
                 if i not in cand]
    texts = [ax.text(lfc[i], lp[i], _short(run.labels.get(i, i), 16), fontsize=7.2,
                     style="italic", color=INK, zorder=6) for i in cand]
    try:
        from adjustText import adjust_text
        obs = s | lfc.index.isin(cand)
        adjust_text(texts, x=lfc[obs].values, y=lp[obs].values, ax=ax,
                    expand=(1.25, 1.6), force_text=(0.4, 0.6), iter_lim=600,
                    ensure_inside_axes=True,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    except Exception:                                                 # noqa: BLE001
        pass
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("−log10(p.adj)" if y_padj else "−log10(p)")
    _title(ax, run.short(c), f"{int(up.sum())} up · {int(dn.sum())} down · "
           + ("colour = functional theme" if groups else "colour = direction")
           + " · names = most marked proteins")
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.11), ncol=2, fontsize=7.5,
              handletextpad=0.3, columnspacing=1.5)
    up_names = [run.labels.get(i, i) for i in score[up & named].sort_values(ascending=False).index[:5]]
    dn_names = [run.labels.get(i, i) for i in score[dn & named].sort_values(ascending=False).index[:5]]
    facts = {"contrast": run.short(c), "up": int(up.sum()), "down": int(dn.sum()),
             "up_names": up_names, "down_names": dn_names,
             "themes": [g[0] for g in groups]}
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", c)[:60]
    return _save(fig, out_dir, f"05_volcano_{idx:02d}_{safe}"), facts


# ---------------------------------------------------------------------------
# 6. Protéines emblématiques
# ---------------------------------------------------------------------------
def fig_spotlight(run, design, prog, out_dir, n=8):
    from .labels import informative
    if run.abundance is None or run.samples is None:
        return None
    lfc = run.de[[c + "_diff" for c in run.contrasts]].abs().max(axis=1)
    picks = []
    if prog is not None:
        for rnd in range(3):
            for p in range(1, prog.k + 1):
                members = prog.labels.index[prog.labels == p]
                cand = [i for i in lfc.loc[members].sort_values(ascending=False).index
                        if informative(run.labels.get(i, "")) and i in run.abundance.index
                        and i not in [x[0] for x in picks]]
                if len(cand) > rnd:
                    picks.append((cand[rnd], p))
                if len(picks) >= n:
                    break
            if len(picks) >= n:
                break
    if len(picks) < n:
        sig = np.zeros(len(run.de), bool)
        for c in run.contrasts:
            sig |= run.sig(c).fillna(False).values
        for i in lfc[sig].sort_values(ascending=False).index:
            if len(picks) >= n:
                break
            if i in run.abundance.index and informative(run.labels.get(i, "")) \
                    and i not in [x[0] for x in picks]:
                picks.append((i, None))
    if not picks:
        return None
    samples = run.samples.set_index("label")
    labs = [l for l in run.abundance.columns]
    cond = samples.loc[labs, "condition"].values
    gcol = group_colors(design)
    lines = design.ordered and len(design.groups or ["all"]) <= 6
    many = (not lines) and len(design.conditions) > 12
    ncol = 2 if many else (4 if len(picks) > 4 else len(picks))
    nrow = int(np.ceil(len(picks) / ncol))
    w = min(7.5, 2.0 + 0.16 * len(design.conditions)) if many else 2.85
    fig, axes = plt.subplots(nrow, ncol, figsize=(w * ncol, (2.3 if many else 2.6) * nrow + 1.0),
                             squeeze=False, sharex=bool(many))
    for ax, (pid, p) in zip(axes.flat, picks):
        v = run.abundance.loc[pid, labs].values.astype(float)
        seen = (run.observed.loc[pid, labs].values == 1) if (
            run.observed is not None and pid in run.observed.index) else np.isfinite(v)
        if lines:
            lv = design.levels
            for g in (design.groups or ["all"]):
                col = gcol.get(g, CAT[0])
                m = np.array([(design.group or {}).get(c, "all") == g for c in cond])
                xs = np.array([lv.index(design.level[c]) for c in cond[m]], float)
                off = (list(design.groups or ["all"]).index(g) - 0.5) * 0.12 if design.group else 0
                ax.scatter(xs[seen[m]] + off, v[m][seen[m]], s=13, color=col, alpha=0.5, lw=0)
                ax.scatter(xs[~seen[m]] + off, v[m][~seen[m]], s=13, facecolors="white",
                           edgecolors=col, lw=0.8)
                mean = pd.Series(v[m]).groupby(xs).mean()
                ax.plot(mean.index, mean.values, color=col, lw=1.8, marker="o", ms=3.5,
                        mec="white", mew=0.7, label=_short(g, 16) if design.group else None)
            ax.set_xticks(range(len(lv)), [design.label(l) for l in lv], fontsize=7.5)
        else:
            conds = design.conditions
            xs = np.array([conds.index(c) for c in cond], float)
            cols = [gcol.get(design.group[c], CAT[0]) if design.group and len(design.groups) <= 3
                    else CAT[0] for c in cond]
            cols = np.array(cols)
            ax.scatter(xs[seen], v[seen], s=13, c=cols[seen], alpha=0.6, lw=0)
            ax.scatter(xs[~seen], v[~seen], s=13, facecolors="white", edgecolors=cols[~seen], lw=0.8)
            mean = pd.Series(v).groupby(xs).mean()
            ax.scatter(mean.index, mean.values, marker="_", s=180, color=INK, lw=1.5, zorder=3)
            ax.set_xticks(range(len(conds)),
                          [_short(design.display(c), 16 if many else 12) for c in conds],
                          rotation=90 if many else 45, ha="center" if many else "right",
                          fontsize=6 if len(conds) > 24 else 7)
        name = str(run.labels.get(pid, pid))
        ax.set_title(_short(name, 20), fontsize=10, style="italic", pad=14)
        desc = ""
        if "First.Protein.Description" in run.de.columns:
            desc = re.sub(r"\s*(OS=|\[).*$", "", str(run.de.at[pid, "First.Protein.Description"]))
        if not desc or desc.lower().startswith(name.rstrip("…").lower()[:12]):
            desc = str(pid)
        ax.text(0, 1.02, _short((f"P{p} · " if p else "") + desc, 40), transform=ax.transAxes,
                fontsize=6.8, color=INK2, va="bottom")
    for ax in axes.flat[len(picks):]:
        ax.set_visible(False)
    if many:   # axe x partagé : noms de conditions sur la dernière ligne visible seulement
        for i, ax in enumerate(axes.flat[:len(picks)]):
            if i + ncol < len(picks):
                ax.tick_params(labelbottom=False)
    for ax in axes[:, 0]:
        ax.set_ylabel(f"log2 abundance ({run.abundance_kind})", fontsize=8)
    if lines and design.group:
        axes.flat[0].legend(fontsize=7, loc="best")
    hollow = " (hollow = not detected, DPC estimate)" if run.observed is not None else ""
    _suptitle(fig, "Spotlight: the proteins behind the programmes",
              f"dots = replicates{hollow} · {'line' if lines else 'dash'} = mean · "
              "most marked named protein of each programme")
    fig.tight_layout(rect=[0, 0, 1, 0.985], h_pad=2.2)
    return _save(fig, out_dir, "06_spotlight"), {"proteins": [run.labels.get(i, i) for i, _ in picks]}
