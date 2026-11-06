"""Bout en bout sur des plans variés : catégoriel, nombreux niveaux, série
temporelle, deux conditions, sortie minimale (une seule feuille)."""
import os
import re

import pytest

from fixture import make_deimos_run
from insight.figures import fig_programmes, style
from insight.factors import detect_design
from insight.io import load_run
from insight.programmes import build_programmes
from insight.run import run_insight
from insight.themes import build_themes


def _files(out):
    return sorted(f for f in os.listdir(out) if f.endswith(".png"))


def _no_external(html):
    return not re.search(r"""(src|href)=["']https?://""", html)


def test_categorical_project(tmp_path, capsys):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat", "Cold"])
    res = run_insight(root)
    assert "[WARN]" not in capsys.readouterr().out
    files = _files(res["out_dir"])
    for p in ("01_landscape", "02_condition_map", "03_programmes", "04_themes", "06_spotlight"):
        assert f"{p}.png" in files, p
    assert any(f.startswith("05_volcano_") for f in files)
    html = open(res["report"], encoding="utf-8").read()
    assert "no order" in html and "Largest step" not in html
    assert _no_external(html)


def test_categorical_layout_is_bars(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat"])
    style()
    run = load_run(root, verbose=False)
    d = detect_design(run.conditions, run.samples)
    prog = build_programmes(run, d, build_themes(run.go))
    assert prog is not None and 3 <= prog.k <= 8
    _, facts = fig_programmes(run, d, prog, str(tmp_path))
    assert facts["layout"] == "bars"


def test_many_conditions_use_heatmap(tmp_path):
    tissues = [f"Tissue{chr(65 + i)}" for i in range(18)]
    root = make_deimos_run(str(tmp_path / "run"), tissues, n_rep=2, n_prot=500)
    style()
    run = load_run(root, verbose=False)
    d = detect_design(run.conditions, run.samples)
    assert not d.ordered
    prog = build_programmes(run, d, build_themes(run.go))
    _, facts = fig_programmes(run, d, prog, str(tmp_path))
    assert facts["layout"] == "heatmap"
    assert os.path.exists(run_insight(root)["report"])


def test_time_course_uses_lines(tmp_path):
    conds = [f"{g}_D{t}" for g in ("WT", "KO") for t in (0, 3, 7, 14)]
    root = make_deimos_run(str(tmp_path / "run"), conds,
                           contrasts=[(f"{g}_D{t}", f"{g}_D0") for g in ("WT", "KO")
                                      for t in (3, 7, 14)] + [("KO_D7", "WT_D7")])
    style()
    run = load_run(root, verbose=False)
    d = detect_design(run.conditions, run.samples)
    assert d.ordered and d.groups == ["WT", "KO"] and d.level_name == "time (d)"
    prog = build_programmes(run, d, build_themes(run.go))
    _, facts = fig_programmes(run, d, prog, str(tmp_path))
    assert facts["layout"] == "lines"
    html = open(run_insight(root)["report"], encoding="utf-8").read()
    assert "Largest step" in html


def test_crossed_design_with_many_groups(tmp_path, capsys):
    """Type Hordeum : 4 génotypes × lumière × temps = 8 groupes ordonnés. Trop de
    groupes pour des trajectoires : centroïdes colorés par niveau, carte de chaleur."""
    gens = ["Tron", "M1", "M2", "M3"]
    conds = [f"{g}_T0_Light" for g in gens] + [f"{g}_T{t}_{l}" for g in gens for t in (3, 7)
                                               for l in ("Light", "Dark")]
    contr = [(f"{g}_T0_Light", "Tron_T0_Light") for g in gens[1:]] + \
        [(f"{g}_T{t}_Light", f"{g}_T{t}_Dark") for g in gens for t in (3, 7)]
    root = make_deimos_run(str(tmp_path / "run"), conds, n_rep=2, n_prot=500, contrasts=contr)
    style()
    run = load_run(root, verbose=False)
    d = detect_design(run.conditions, run.samples)
    assert d.ordered and len(d.groups) == 8 and d.level_name == "time"
    prog = build_programmes(run, d, build_themes(run.go))
    _, facts = fig_programmes(run, d, prog, str(tmp_path))
    assert facts["layout"] == "heatmap"
    capsys.readouterr()
    html = open(run_insight(root)["report"], encoding="utf-8").read()
    assert "[WARN]" not in capsys.readouterr().out
    assert "Largest step" not in html


def test_two_conditions(tmp_path, capsys):
    root = make_deimos_run(str(tmp_path / "run"), ["Grp_Male4", "Grp_Female4"], n_rep=4)
    res = run_insight(root)
    out = capsys.readouterr().out
    assert "[WARN]" not in out and "no order" in out
    files = _files(res["out_dir"])
    assert "03_programmes.png" not in files          # < 3 conditions : pas de programmes
    assert "01_landscape.png" in files and "06_spotlight.png" in files


def test_minimal_workbook(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), ["A", "B", "C"], go=False, matrices=False,
                           pca=False, manifest=False)
    res = run_insight(root)
    files = _files(res["out_dir"])
    assert "01_landscape.png" in files
    assert any(f.startswith("05_volcano_") for f in files)
    assert "04_themes.png" not in files and "02_condition_map.png" not in files
    assert os.path.exists(res["report"])


def test_themes_merge_redundant_terms(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat"])
    run = load_run(root, verbose=False)
    th = build_themes(run.go)
    # deux termes par programme simulé, qui partagent ~85 % de leurs protéines
    assert len(th.names) < th.terms.term.nunique()
    assert th.terms.groupby("theme").term.nunique().max() >= 2


def test_selection_follows_deimos(tmp_path):
    """Protéines = celles de Deimos : clés structurées > texte > last_config > défauts."""
    import numpy as np
    import pandas as pd
    from insight.io import _rule
    man = pd.DataFrame([("stats", "volcano_rule", "p.adj < 0.05 & ratio >= 1.5"),
                        ("stats", "volcano_p_column", "p.val"),
                        ("stats", "volcano_p_thresh", "0.01"),
                        ("stats", "volcano_ratio_min", "2.0")],
                       columns=["section", "key", "value"])
    assert _rule(man, {}) == (False, 0.01, 2.0, "Deimos manifest")
    assert _rule(man.iloc[:1], {})[:3] == (True, 0.05, 1.5)
    assert _rule(None, {"volcano_use_padj": True, "volcano_p_thresh": 0.1,
                        "volcano_ratio_min": 2})[:3] == (True, 0.1, 2.0)
    assert _rule(None, {})[:3] == (False, 0.05, 1.5)            # défauts de Deimos
    root = make_deimos_run(str(tmp_path / "run"), ["A", "B", "C"])
    run = load_run(root, verbose=False)
    assert run.rule_source == "Deimos manifest" and run.use_padj and run.ratio_min == 2.0
    for c in run.contrasts:                                       # = volcano_status
        p, d = run.de[c + "_p.adj"].to_numpy(), run.de[c + "_diff"].to_numpy()
        assert int(run.sig(c).sum()) == int(((p < 0.05) & (np.abs(d) >= 1.0)).sum())


def test_go_terms_kept_at_q_005(tmp_path):
    import pandas as pd
    root = make_deimos_run(str(tmp_path / "run"), ["A", "B", "C"])
    xl = os.path.join(root, "ProteomicAnalysis_Results_test.xlsx")
    sheets = pd.read_excel(xl, sheet_name=None)
    go = [n for n in sheets if n.startswith("GO_")][0]
    t = sheets[go].copy()
    t.loc[0, "p_value"] = 0.2                                     # terme non significatif
    sheets[go] = t
    with pd.ExcelWriter(xl) as w:
        for n, d in sheets.items():
            d.to_excel(w, sheet_name=n, index=False)
    run = load_run(root, verbose=False)
    allp = pd.concat(run.go.values()).p_value
    assert (allp < 0.05).all() and len(allp) == sum(
        (pd.to_numeric(sheets[n].p_value) < 0.05).sum() for n in sheets if n.startswith("GO_"))
