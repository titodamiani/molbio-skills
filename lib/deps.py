"""Check the packages a skill needs, and say how to install them."""
import importlib
import sys
from pathlib import Path

REQUIREMENTS = Path(__file__).resolve().parents[1] / "requirements.txt"

# Import name -> pip name, for the ones that differ.
PIP_NAMES = {"Bio": "biopython", "primer3": "primer3-py"}


def require(*modules):
    """Stop with the exact pip command when something is missing.

    Nothing is installed here. A skill that installs packages by itself makes
    the pinned versions in requirements.txt meaningless, because pip fetches
    the newest release instead of the pin.
    """
    missing = [PIP_NAMES.get(module, module) for module in modules
               if not _importable(module)]
    if missing:
        sys.exit(f"missing: {', '.join(missing)}\n"
                 f"install with:\n"
                 f"    python3 -m pip install -r {REQUIREMENTS}")


def _importable(module):
    try:
        importlib.import_module(module)
    except ImportError:
        return False
    return True
