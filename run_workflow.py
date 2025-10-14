import argparse
import subprocess
import json
import tempfile
import os

def main():
    parser = argparse.ArgumentParser(description="Wrapper for Snakemake workflow")

    # Workflow config
    parser.add_argument("--directories", nargs="+", required=True, help="List of input directories")
    parser.add_argument("--metabolome", action="store_true", help="Enable metabolome analysis")
    parser.add_argument("--genome", action="store_true", help="Enable genome analysis")
    parser.add_argument("--e-weights", action="store_true", help="Enable E weights")
    parser.add_argument("--n-weights", action="store_true", help="Enable N weights")
    parser.add_argument("--include-orgs", action="store_true", help="Include organisms")
    parser.add_argument("--abundance-col", type=str, default="Abundance_RPKs", help="Column name for abundance")
    parser.add_argument('--uri', type=str, required=True, help='Neo4j URI instance')
    parser.add_argument('--user', type=str, required=True, help='Neo4j username for instance')
    parser.add_argument('--p', type=str, required=True, help='Neo4j password for instance')

    # Snakemake execution options
    parser.add_argument("--cores", type=int, default=1, help="Number of cores to use")
    parser.add_argument("--profile", type=str, help="Snakemake profile to use")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Perform a dry-run")

    args = parser.parse_args()

    # Convert directories list to comma-separated string
    directories_str = ",".join(args.directories)

    # Prepare config dict
    config_args = {
        "directories": directories_str,
        "metabolome": args.metabolome,
        "genome": args.genome,
        "e_weights": args.e_weights,
        "n_weights": args.n_weights,
        "include_orgs": args.include_orgs,
        "abundance_col": args.abundance_col,
        "uri": args.uri,
        "user": args.user,
        "p": args.p,
    }

    # ✅ Write config to a temporary JSON file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tmp:
        json.dump(config_args, tmp, indent=2)
        tmp_path = tmp.name

    # Base Snakemake command using --configfile (safe for all characters)
    cmd = f"snakemake --snakefile Snakefile --cores {args.cores} --configfile {tmp_path}"

    # Optional flags
    if args.profile:
        cmd += f" --profile {args.profile}"
    if args.dry_run:
        cmd += " -n"

    print("🚀 Running Snakemake command:")
    print(cmd)
    print(f"🧾 Using config file: {tmp_path}")

    try:
        subprocess.run(cmd, shell=True, check=True)
    finally:
        # Clean up temp config file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

if __name__ == "__main__":
    main()