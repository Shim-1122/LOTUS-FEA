import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def febio_available():
    return bool(os.environ.get("RLFEA_FEBIO") or shutil.which("febio4"))


requires_febio = pytest.mark.skipif(not febio_available(), reason="FEBio executable not available")
