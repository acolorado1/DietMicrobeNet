from pathlib import Path
import subprocess
import sys
from .downloader import download_data


ROOT = Path(__file__).resolve().parents[3]


def run_GetFoods():

    app = ROOT / "src" / "get_foods.py"

    subprocess.run(
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            str(app)
        ],
        check=True
    )


def run_workflow():

    workflow = ROOT / "run_workflow.py"

    subprocess.run(
        [sys.executable, str(workflow)],
        check=True
    )


def run_DM_GraphComparison():

    graph = ROOT / "src" / "GraphComparison.py"

    subprocess.run(
        [sys.executable, str(graph)],
        check=True
    )


def run_DMHost_GraphComparison():

    graph = (
        ROOT
        / "src"
        / "Host"
        / "host_GraphComparison.py"
    )

    subprocess.run(
        [sys.executable, str(graph)],
        check=True
    )