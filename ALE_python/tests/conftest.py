import pytest
import os
import subprocess
import tempfile
import shutil

@pytest.fixture
def project_root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

@pytest.fixture
def example_data(project_root):
    return os.path.join(project_root, "example_data")

@pytest.fixture
def cpp_bin(project_root):
    return os.path.join(project_root, "build", "bin")

@pytest.fixture
def species_tree(example_data):
    return os.path.join(example_data, "S.tree")

@pytest.fixture
def gene_trees(example_data):
    return os.path.join(example_data, "HBG745965_real.1.treelist")

@pytest.fixture
def tmp_dir():
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)
