import json
from pathlib import Path
import pytest

P = Path(__file__).resolve().parents[1] / "results" / "structural_similarity.json"


@pytest.mark.skipif(not P.exists(), reason="result not generated")
def test_structural_similarity_result():
    j = json.loads(P.read_text())
    assert j["n_pairs"] > 0 and 0 < j["synergy_base_rate"] < 1
    assert set(j["synergy_rate_by_quartile"]) == {"1", "2", "3", "4"}
    a = j["adjusted_logit"]
    assert 0 <= a["p"]["tanimoto"] <= 1 and a["n_perm"] == 200


def test_tanimoto_identity():
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    g = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    f = g.GetFingerprint(Chem.MolFromSmiles("CCO"))
    assert DataStructs.TanimotoSimilarity(f, f) == 1.0
