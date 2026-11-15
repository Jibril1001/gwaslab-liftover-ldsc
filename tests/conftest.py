import sys
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import synthetic  # noqa: E402

T2T_CHAIN = ROOT / "data" / "reference" / "grch38-chm13v2.chain"


@pytest.fixture(scope="session")
def manifest(tmp_path_factory):
    """Synthetic BBJ-style sumstats + LDSC references, generated once per session."""
    return synthetic.generate(tmp_path_factory.mktemp("synthetic"))


@pytest.fixture(scope="session")
def truth(manifest):
    return manifest["truth"]
