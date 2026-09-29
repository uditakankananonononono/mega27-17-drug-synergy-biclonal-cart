"""Hermetic tests for cart-target-rank on tiny synthetic fixtures."""
import json, pathlib, subprocess, sys, zipfile, os

REPO = pathlib.Path(__file__).resolve().parents[1]
ENV = {**os.environ, "PYTHONPATH": str(REPO / "src")}


def make_data(tmp):
    d = tmp / "data"; d.mkdir()
    (d / "rna_tissue_consensus.tsv").write_text(
        "Gene\tGene name\tTissue\tnTPM\n"
        "ENSG1\tGENEA\tlung\t10\nENSG1\tGENEA\theart muscle\t20\n"
        "ENSG2\tGENEB\tlung\t5\nENSG2\tGENEB\tskin\t50\n")
    (d / "subcellular_location.tsv").write_text(
        "Gene\tMain location\tAdditional location\nENSG1\tPlasma membrane\t\nENSG2\tNucleoplasm\t\n")
    (d / "uniprot_cellmembrane.tsv").write_text("Entry\tGene Names\nP1\tGENEB\n")
    rows = ["Gene\tGene name\tCancer\tFPKM"]
    for i, v in enumerate([5, 15, 25, 100]):
        rows.append(f"ENSG1\tGENEA\tBRCA\t{v}")
    for v in [1, 2]:
        rows.append(f"ENSG2\tGENEB\tBRCA\t{v}")
    with zipfile.ZipFile(d / "rna_cancer_sample.tsv.zip", "w") as zf:
        zf.writestr("rna_cancer_sample.tsv", "\n".join(rows) + "\n")
    return d


def test_rank_genes(tmp_path):
    sys.path.insert(0, str(REPO / "src"))
    from cart.cli import rank_genes
    d = make_data(tmp_path)
    rows = rank_genes(["GENEA", "GENEB"], "BRCA", str(d))
    a = rows["GENEA"]
    assert a["surface"] is True and a["n_samples"] == 4
    assert a["normal_max_nTPM"] == 20.0          # vital tissues only
    assert a["fanm"] == 0.5                       # 2 of 4 samples above 20
    b = rows["GENEB"]
    assert b["surface"] is True                   # via UniProt union
    assert b["normal_max_nTPM"] == 5.0            # skin excluded (non-vital)
    assert b["fanm"] == 0.0


def test_cli_json(tmp_path):
    d = make_data(tmp_path)
    out = tmp_path / "out.json"
    r = subprocess.run([sys.executable, "-m", "cart.cli", "GENEA", "--cancer", "BRCA",
                        "--data-dir", str(d), "--json", str(out)],
                       capture_output=True, text=True, env=ENV, timeout=120)
    assert r.returncode == 0, r.stderr
    res = json.loads(out.read_text())
    assert res["genes"]["GENEA"]["fanm"] == 0.5


def test_report_mode(tmp_path):
    """Hermetic --report test: synthetic corpus join + boundary note for out-of-corpus genes."""
    sys.path.insert(0, str(REPO / "src"))
    from cart.cli import main, load_report_corpus, BOUNDARY_NOTE
    d = make_data(tmp_path)
    res = tmp_path / "results"; res.mkdir()
    (res / "pharos_per_gene.csv").write_text(
        "gene,group,status,tdl,novelty,n_drugs,n_ligands,publication_count,generif_count,gwas_total,ppi_total\n"
        "GENEA,gated,ok,Tclin,0.5,2,3,10,5,1,7\n")
    out = tmp_path / "out.json"
    main(["GENEA", "GENEB", "--cancer", "BRCA", "--data-dir", str(d),
          "--report", "--results-dir", str(res), "--json", str(out)])
    d2 = json.load(open(out))
    assert d2["genes"]["GENEA"]["audit"]["pharos"]["tdl"] == "Tclin"
    assert d2["genes"]["GENEA"]["corpus"] == "188-gene study corpus"
    assert d2["genes"]["GENEB"]["audit"] == {}
    assert d2["genes"]["GENEB"]["corpus"] == BOUNDARY_NOTE
    corpus = load_report_corpus(res)
    assert corpus["GENEA"]["pharos"]["n_drugs"] == "2"
