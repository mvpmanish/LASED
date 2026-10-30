"""The README quick start and the example scripts must keep working."""
import pathlib
import re
import runpy

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_readme_quickstart():
    readme = (ROOT/"README.md").read_text()
    block = re.search(r"## Quick start\s+```python\n(.*?)```", readme, re.S).group(1)
    namespace = {}
    exec(block, namespace)
    assert namespace["rho"].shape == (501, 4, 4)
    rho = namespace["rho"]
    assert np.abs(rho - rho.conj().transpose(0, 2, 1)).max() < 1e-12


@pytest.mark.parametrize("script", sorted((ROOT/"examples").glob("*.py")), ids=lambda p: p.name)
def test_example_runs(script, capsys):
    runpy.run_path(str(script), run_name="__main__")
    assert capsys.readouterr().out
