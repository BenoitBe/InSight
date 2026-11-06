"""Export des tables (insight_tables.xlsx, insight_programmes.tsv) : feuilles, cohérence
avec les programmes et thèmes des figures, part observée (DPC), plans sans programme."""
import os

import numpy as np
import pandas as pd

from fixture import make_deimos_run
from insight.run import run_insight
from insight.tables import TSV_NAME, XLSX_NAME


def _add_dpc(root, seed=1):
    """Ajoute dpc_abundance.tsv / dpc_observed.tsv : 20 protéines jamais vues dans une
    condition, pour vérifier min_observed_fraction."""
    md = os.path.join(root, "matrices")
    A = pd.read_csv(os.path.join(md, "log2_imputed.tsv"), sep="\t")
    s = pd.read_csv(os.path.join(md, "samples.tsv"), sep="\t")
    A.to_csv(os.path.join(md, "dpc_abundance.tsv"), sep="\t", index=False)
    O = A.copy()
    labs = list(s.label)
    O[labs] = (np.random.default_rng(seed).random((len(O), len(labs))) > 0.05).astype(int)
    first = s.condition.iloc[0]
    O.loc[:19, list(s.label[s.condition == first])] = 0
    O.to_csv(os.path.join(md, "dpc_observed.tsv"), sep="\t", index=False)
    return list(O.protein[:20]), first


def test_tables_content(tmp_path, capsys):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat", "Cold"])
    gone, first = _add_dpc(root)
    res = run_insight(root)
    assert "[WARN]" not in capsys.readouterr().out
    assert res["tables"] and os.path.basename(res["tables"]) == XLSX_NAME
    x = pd.ExcelFile(res["tables"])
    assert x.sheet_names == ["README", "Facts", "Programmes", "Programme_summary", "Themes",
                             "Theme_terms", "Theme_x_contrast"]
    pr = pd.read_excel(x, "Programmes")
    summ = pd.read_excel(x, "Programme_summary")
    # chaque protéine une fois ; tailles cohérentes avec le résumé
    assert pr.protein.is_unique
    assert pr.programme.value_counts().to_dict() == dict(zip(summ.programme, summ.n_proteins))
    assert pr.r_to_programme.between(-1, 1).all() and pr.silhouette.between(-1, 1).all()
    assert (pr.n_significant >= 1).all()                  # protéines de la sélection Deimos
    for c in ["Ctrl", "DrugA", "DrugB", "Heat", "Cold"]:
        assert f"log2_{c}" in pr.columns and f"z_{c}" in pr.columns
    z = pr[[c for c in pr.columns if c.startswith("z_")]]
    assert np.allclose(z.mean(axis=1), 0, atol=1e-2)
    # part observée : les protéines jamais vues dans la 1re condition
    sel = pr[pr.protein.isin(gone)]
    assert len(sel) and (sel.min_observed_fraction == 0).all()
    assert (sel.least_observed_in == first).all()
    # thèmes : mêmes identifiants dans les trois feuilles
    th, tt, tx = (pd.read_excel(x, s) for s in ("Themes", "Theme_terms", "Theme_x_contrast"))
    assert set(tt.theme_id) == set(th.theme_id) == set(tx.theme_id)
    assert th.best_fdr.is_monotonic_increasing and (tt.best_fdr < 0.05).all()
    assert (th.n_terms == tt.groupby("theme_id").size().reindex(th.theme_id).values).all()
    facts = pd.read_excel(x, "Facts")
    assert "Response programmes" in set(facts.section)
    assert not facts.fact.str.contains("<").any()          # texte sans balises HTML
    readme = pd.read_excel(x, "README", header=None)
    assert readme.iloc[:, 1].astype(str).str.contains("p.adj < 0.05").any()
    tsv = pd.read_csv(os.path.join(res["out_dir"], TSV_NAME), sep="\t")
    assert list(tsv.protein) == list(pr.protein)
    html = open(res["report"], encoding="utf-8").read()
    assert XLSX_NAME in html


def test_tables_without_programmes(tmp_path):
    """Deux conditions : pas de programme ni de TSV, mais Facts et thèmes."""
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "Treated"])
    res = run_insight(root)
    x = pd.ExcelFile(res["tables"])
    assert "Programmes" not in x.sheet_names and "Facts" in x.sheet_names
    assert not os.path.exists(os.path.join(res["out_dir"], TSV_NAME))


def test_tables_minimal_workbook(tmp_path):
    """Classeur réduit à Differential_Expression : tables quand même (Facts, README)."""
    root = make_deimos_run(str(tmp_path / "run"), ["A", "B", "C"], go=False, matrices=False,
                           pca=False, manifest=False)
    res = run_insight(root)
    x = pd.ExcelFile(res["tables"])
    assert x.sheet_names[:2] == ["README", "Facts"]
