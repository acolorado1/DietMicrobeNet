from pathlib import Path
import requests


FILES = {
    "CompoundExternalDescriptor.csv": "https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/ESXx7vpypQFOt4iVv6x-ErkBykpAVS1fppQjYZkrxkDnAA?download=1",
    "Content.csv": "https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EYJUYQWmY9VDlYZIAXpzpvEBzhrnViFZQjrikXIla_aPPg?download=1",
    "Food.csv": "https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EXyRAlYs1htNlcwz5T67BxQBGO7HfOjmfIBlkOydM0BIAw?download=1",
    "hmdb.csv": "https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EbY2fD3JTcNLomKFqQhY5jABAXN-60A80PmkngRynazocg?download=1",
    "AllFood/food_meta.csv": "https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EZ1pyHd616RFkR9zG6kenuoBhZDroHYTbaGmEfwpxFOHLg?download=1",
}


def download_data():

    root = Path.cwd()

    data_dir = root / "Data"
    data_dir.mkdir(exist_ok=True)

    for filename, url in FILES.items():

        outfile = data_dir / filename

        outfile.parent.mkdir(
            parents=True,
            exist_ok=True
        )

        if outfile.exists():
            print(f"Skipping {filename}")
            continue

        print(f"Downloading {filename}")

        r = requests.get(url)
        r.raise_for_status()

        with open(outfile, "wb") as f:
            f.write(r.content)