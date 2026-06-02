from pathlib import Path
import subprocess
import sys

from .downloader import download_data  # noqa: F401  (re-exported for entry point)

# ---------------------------------------------------------------------------
# Project-root resolution
# ---------------------------------------------------------------------------

_SENTINEL_FILES = ("run_workflow.py", "Snakefile", "pyproject.toml")


def _find_root() -> Path:
    """Walk up from the current working directory looking for the project root.

    The root is identified by the presence of known top-level files.  This
    approach works whether the package was installed in editable mode (*pip
    install -e .*) or as a regular wheel, because it searches from wherever
    the user invoked the command rather than relying on the package's
    installation path.

    Raises
    ------
    FileNotFoundError
        When no ancestor directory contains the expected sentinel files.
    """
    for candidate in [Path.cwd(), *Path.cwd().parents]:
        if all((candidate / s).exists() for s in _SENTINEL_FILES):
            return candidate

    raise FileNotFoundError(
        "Could not locate the DietMicrobeNet project root.\n"
        "Please run this command from inside the cloned repository "
        "(the directory that contains run_workflow.py and pyproject.toml)."
    )


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def run_GetFoods() -> None:
    """Launch the FooDB food-selection Streamlit app."""
    root = _find_root()
    app = root / "src" / "get_foods.py"
    _require_file(app, "get_foods.py")

    result = subprocess.run(
        [sys.executable, "-m", "streamlit", "run", str(app)] + sys.argv[1:],
    )
    sys.exit(result.returncode)


def run_workflow() -> None:
    """Run the main DietMicrobeNet Snakemake workflow."""
    root = _find_root()
    workflow = root / "run_workflow.py"
    _require_file(workflow, "run_workflow.py")

    result = subprocess.run(
        [sys.executable, str(workflow)] + sys.argv[1:],
    )
    sys.exit(result.returncode)


def run_DM_GraphComparison() -> None:
    """Run DietMicrobe graph comparison script"""
    root = _find_root()
    graph = root / "src" / "GraphComparison.py"
    _require_file(graph, "GraphComparison.py")

    result = subprocess.run(
        [sys.executable, str(graph)] + sys.argv[1:],
    )
    sys.exit(result.returncode)


def run_DMHost_GraphComparison() -> None:
    """Run DietMicrobeHost graph comparison script """
    root = _find_root()
    graph = root / "src" / "Host" / "host_GraphComparison.py"
    _require_file(graph, "host_GraphComparison.py")

    result = subprocess.run(
        [sys.executable, str(graph)] + sys.argv[1:],
    )
    sys.exit(result.returncode)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _require_file(path: Path, label: str) -> None:
    """Raise a clear error when an expected script is missing."""
    if not path.exists():
        raise FileNotFoundError(
            f"Expected script not found: {path}\n"
            f"Make sure '{label}' is present in the repository and that you "
            "are running the command from the project root."
        )