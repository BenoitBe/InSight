"""Orchestration : lit un run Deimos, construit thèmes et programmes, produit les
figures et la page récit. Chaque étape est protégée : une figure qui échoue est
signalée et n'empêche pas les autres."""
from __future__ import annotations

import html
import os
import traceback

from . import __version__
from .factors import detect_design
from .figures import (fig_landscape, fig_map, fig_programmes, fig_spotlight, fig_themes,
                      fig_volcano, style)
from .io import load_run
from .programmes import build_programmes
from .report import write_report
from .tables import write_tables
from .themes import build_themes


def _safe(label, fn, *a, **kw):
    try:
        return fn(*a, **kw)
    except Exception as e:                                            # noqa: BLE001
        print(f"  [WARN] InSight: {label} skipped ({type(e).__name__}: {e})")
        if os.environ.get("INSIGHT_DEBUG"):
            traceback.print_exc()
        return None


def _names(lst, n=5):
    return ", ".join(f"<em>{html.escape(str(x))}</em>" for x in lst[:n]) or "—"


def run_insight(deimos_out: str, out_dir: str | None = None, xlsx: str | None = None,
                reference_fasta: str | None = None, k: int | None = None,
                volcanos: list | None = None, max_volcanos: int = 4,
                design_overrides: dict | None = None, title: str | None = None) -> dict:
    style()
    out_dir = out_dir or os.path.join(deimos_out, "insight")
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n[InSight {__version__}] Biological reading of {deimos_out}")
    run = load_run(deimos_out, xlsx=xlsx, reference_fasta=reference_fasta)
    design = detect_design(run.conditions, run.samples, design_overrides)
    kind = (f"{len(design.groups) or 1} group(s) × ordered {design.level_name} "
            f"({', '.join(design.label(l) for l in design.levels)})" if design.ordered
            else f"{len(design.conditions)} conditions, no order")
    print(f"  [InSight] design: {kind} [{design.source}]")
    themes = _safe("themes", build_themes, run.go)
    prog = _safe("programmes", build_programmes, run, design, themes, k)

    sections, written = [], []

    r = _safe("landscape", fig_landscape, run, out_dir)
    if r:
        paths, f = r; written += paths
        b = [f"{f['n_contrasts']} comparisons; {f['empty']} without any significant protein."]
        b += [f"<b>{html.escape(n)}</b>: {u} up, {d} down." for n, u, d in f["top"]]
        sections.append(dict(id="landscape", title="How much changes", png=paths[0], bullets=b))

    r = _safe("condition map", fig_map, run, design, out_dir)
    if r:
        paths, f = r; written += paths
        b = []
        v1, v2 = f.get("var", (None, None))
        if v1 is not None:
            b.append(f"The first two components carry {v1:.0f} % and {v2:.0f} % of the variance.")
        if f.get("farthest"):
            a, c = f["farthest"]
            b.append(f"Most distant conditions: <b>{html.escape(a)}</b> and <b>{html.escape(c)}</b>.")
        for g, a, c, d in f.get("steps", [])[:2]:
            b.append(f"Largest step{'' if g == 'all' else ' in ' + html.escape(str(g))}: "
                     f"{html.escape(a)} → {html.escape(c)}.")
        sections.append(dict(id="map", title="Where the conditions sit", png=paths[0], bullets=b))

    r = _safe("programmes", fig_programmes, run, design, prog, out_dir)
    if r:
        paths, _ = r; written += paths
        b = []
        for p in range(1, prog.k + 1):
            n = int((prog.labels == p).sum())
            th = prog.themes.get(p, [])
            ths = "; ".join(html.escape(t[0]) for t in th) if th else "no enriched theme"
            b.append(f"<b>P{p}</b> ({n} proteins): highest in {html.escape(prog.high[p])}, "
                     f"lowest in {html.escape(prog.low[p])} · themes: {ths} · "
                     f"e.g. {_names(prog.names.get(p, []), 4)}.")
        sections.append(dict(id="programmes", title="Response programmes", png=paths[0], bullets=b))

    r = _safe("themes", fig_themes, run, themes, out_dir)
    if r:
        paths, f = r; written += paths
        b = [f"{f['n_terms']} enriched terms grouped into {f['n_themes']} themes "
             f"(terms sharing proteins)."]
        b += [f"<b>{html.escape(n)}</b>: strongest in {html.escape(', '.join(cs))}."
              for n, cs in f["top"]]
        sections.append(dict(id="themes", title="Functional themes", png=paths[0], bullets=b))

    if volcanos:
        vc = [c for c in volcanos if c in run.contrasts]
    else:
        counts = {c: int(run.sig(c).fillna(False).sum()) for c in run.contrasts}
        vc = [c for c, n in sorted(counts.items(), key=lambda x: -x[1]) if n > 0][:max_volcanos]
    for i, c in enumerate(vc, 1):
        r = _safe(f"volcano {c}", fig_volcano, run, c, themes, out_dir, i)
        if r:
            paths, f = r; written += paths
            b = [f"{f['up']} up, {f['down']} down."]
            if f["up_names"]:
                b.append(f"Most marked up: {_names(f['up_names'])}.")
            if f["down_names"]:
                b.append(f"Most marked down: {_names(f['down_names'])}.")
            if f["themes"]:
                b.append("Coloured themes: " + "; ".join(html.escape(t) for t in f["themes"]) + ".")
            sections.append(dict(id=f"volcano{i}", title=f["contrast"], png=paths[0], bullets=b))

    r = _safe("spotlight", fig_spotlight, run, design, prog, out_dir)
    if r:
        paths, f = r; written += paths
        sections.append(dict(id="spotlight", title="Spotlight proteins", png=paths[0],
                             bullets=[f"Shown: {_names(f['proteins'], 8)}."]))

    n_sig = int(sum(run.sig(c).fillna(False) for c in run.contrasts).gt(0).sum()) \
        if run.contrasts else 0
    kpis = [("proteins quantified", len(run.de)), ("comparisons", len(run.contrasts)),
            ("differential (any comparison)", n_sig)]
    if prog is not None:
        kpis.append(("response programmes", prog.k))
    if themes is not None:
        kpis.append(("functional themes", len(themes.names)))
    from .io import GO_Q_THR
    meta = [f"Source: {os.path.basename(run.xlsx)} · design: {kind}",
            f"Proteins: {run.rule_text()} ({run.rule_source}: the Deimos selection) · "
            f"GO terms: q < {GO_Q_THR:g} · abundances: {run.abundance_kind or 'n/a'} · "
            f"protein names: {run.label_source}"]
    tables = _safe("tables", write_tables, out_dir, run, design, prog, themes, sections, kpis,
                   meta, k)
    if tables:
        written.append(tables)
        meta.append(f"Tables: {os.path.basename(tables)} (same folder) — programmes, themes, "
                    f"facts" + (" · programmes also in insight_programmes.tsv" if prog else ""))
    rep = write_report(out_dir, title or "InSight — biological reading", kpis, sections, meta)
    print(f"  [OK] InSight: {len(sections)} figures + report -> {rep}")
    if tables:
        print(f"  [OK] InSight: tables -> {tables}")
    return {"report": rep, "files": written, "out_dir": out_dir, "tables": tables}
