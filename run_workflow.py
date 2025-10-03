import argparse
import subprocess
import shlex

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

    # Snakemake execution options
    parser.add_argument("--snakefile", type=str, default="Snakefile", help="Path to Snakefile")
    parser.add_argument("--cores", type=int, default=1, help="Number of cores to use")
    parser.add_argument("--profile", type=str, help="Snakemake profile to use")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Perform a dry-run")

    args = parser.parse_args()

    # Build config dictionary for Snakemake
    config_args = {
        "directories": args.directories,
        "metabolome": args.metabolome,
        "genome": args.genome,
        "e_weights": args.e_weights,
        "n_weights": args.n_weights,
        "include_orgs": args.include_orgs,
        "abundance_col": args.abundance_col,
    }

    # Convert config dict into CLI string
    config_str = " ".join(
        f'{key}="{value}"' if not isinstance(value, list) else f'{key}="{value}"'
        for key, value in config_args.items()
    )

    # Base snakemake command
    cmd = f"snakemake --snakefile {args.snakefile} --cores {args.cores} --config {config_str}"

    # Optional flags
    if args.profile:
        cmd += f" --profile {args.profile}"
    if args.dry_run:
        cmd += " -n"

    print("Running:", cmd)

    # Run snakemake
    subprocess.run(shlex.split(cmd), check=True)

if __name__ == "__main__":
    main()
