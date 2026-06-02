from pathlib import Path
import sys

import requests

# ---------------------------------------------------------------------------
# Files to download
# ---------------------------------------------------------------------------

FILES: dict[str, str] = {
    "CompoundExternalDescriptor.csv": (
        "https://olucdenver-my.sharepoint.com/:x:/g/personal/"
        "angelasofia_burkhartcolorado_cuanschutz_edu/"
        "ESXx7vpypQFOt4iVv6x-ErkBykpAVS1fppQjYZkrxkDnAA?download=1"
    ),
    "Content.csv": (
        "https://olucdenver-my.sharepoint.com/:x:/g/personal/"
        "angelasofia_burkhartcolorado_cuanschutz_edu/"
        "EYJUYQWmY9VDlYZIAXpzpvEBzhrnViFZQjrikXIla_aPPg?download=1"
    ),
    "Food.csv": (
        "https://olucdenver-my.sharepoint.com/:x:/g/personal/"
        "angelasofia_burkhartcolorado_cuanschutz_edu/"
        "EXyRAlYs1htNlcwz5T67BxQBGO7HfOjmfIBlkOydM0BIAw?download=1"
    ),
    "hmdb.csv": (
        "https://olucdenver-my.sharepoint.com/:x:/g/personal/"
        "angelasofia_burkhartcolorado_cuanschutz_edu/"
        "EbY2fD3JTcNLomKFqQhY5jABAXN-60A80PmkngRynazocg?download=1"
    ),
    "AllFood/food_meta.csv": (
        "https://olucdenver-my.sharepoint.com/:x:/g/personal/"
        "angelasofia_burkhartcolorado_cuanschutz_edu/"
        "EZ1pyHd616RFkR9zG6kenuoBhZDroHYTbaGmEfwpxFOHLg?download=1"
    ),
}

# How long to wait for a server response before giving up (seconds).
_TIMEOUT = 60


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def download_data() -> None:
    """Download required data files into *<cwd>/Data/*.

    Files that already exist on disk are skipped.  Failed downloads are
    reported to stderr but do not abort the remaining downloads.

    Notes
    -----
    Files are written relative to the current working directory so that they
    land in the correct ``Data/`` folder regardless of how the package was
    installed.  Run this command from the repository root (the directory that
    contains ``pyproject.toml``) to ensure the files are placed correctly.
    """
    root = Path.cwd()

    # Warn when the user is likely in the wrong directory.
    if not (root / "pyproject.toml").exists():
        print(
            "Warning: 'pyproject.toml' not found in the current directory.\n"
            "         Data files will be saved to ./Data/ here instead of the\n"
            "         project root.  Consider running from the repository root.",
            file=sys.stderr,
        )

    data_dir = root / "Data"
    data_dir.mkdir(exist_ok=True)

    failed: list[str] = []

    for filename, url in FILES.items():
        outfile = data_dir / filename
        outfile.parent.mkdir(parents=True, exist_ok=True)

        if outfile.exists():
            print(f"  skipping  {filename}  (already exists)")
            continue

        print(f"  downloading  {filename} … ", end="", flush=True)

        try:
            response = requests.get(url, timeout=_TIMEOUT)
            response.raise_for_status()
        except requests.Timeout:
            print("FAILED")
            print(
                f"    Timeout after {_TIMEOUT} s — server did not respond.",
                file=sys.stderr,
            )
            failed.append(filename)
            continue
        except requests.HTTPError as exc:
            print("FAILED")
            print(f"    HTTP error: {exc}", file=sys.stderr)
            failed.append(filename)
            continue
        except requests.RequestException as exc:
            print("FAILED")
            print(f"    Network error: {exc}", file=sys.stderr)
            failed.append(filename)
            continue

        outfile.write_bytes(response.content)
        print("done")

    # Summary
    total = len(FILES)
    succeeded = total - len(failed)
    print(f"\n{succeeded}/{total} file(s) downloaded successfully.")

    if failed:
        print("Failed downloads:", file=sys.stderr)
        for name in failed:
            print(f"  • {name}", file=sys.stderr)
        sys.exit(1)