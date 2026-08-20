import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Generate HTML report of microbe originating compounds.")
parser.add_argument("--patterns", default="default.csv", help="Path to the graph_results.csv")
parser.add_argument("--metabolome", default="default.csv", help="Path to CSV file containing one column of KEGG compounds.")
parser.add_argument("--output", default="filepath.html", help="Path to HTML report file")
args = parser.parse_args()

# Locate the Rmd file relative to this script's own location.
this_dir = os.path.dirname(os.path.abspath(__file__))
rmd_path = os.path.join(this_dir, "MetabolomeComparison_Report.Rmd")

# IMPORTANT: rmarkdown::render() knits with the working directory set to
# the .Rmd file's own directory (this_dir), NOT the directory the user ran
# this script from. Since the .Rmd reads params$patterns / params$metabolome
# using relative-path semantics, we must resolve the user's inputs to
# absolute paths *before* passing them to R - anchored to the user's
# current working directory (where they actually ran this script), not to
# this_dir. Otherwise R will look for these files in the wrong place.
patterns_path = os.path.abspath(args.patterns)
metabolome_path = os.path.abspath(args.metabolome)

# NOTE: unlike the other report scripts, this one treats --output as a
# DIRECTORY and appends a fixed filename, per the original script's
# behavior. Resolved to an absolute path against the user's cwd for the
# same reason as above.
output_dir = os.path.abspath(args.output)
output_file_path = os.path.join(output_dir, "MetabolomeComparison_Report.html")

# Construct the R command as a list (not a shell string) so paths
# containing spaces, quotes, or other special characters can't break
# the command or introduce shell-injection issues.
r_command = [
    "Rscript",
    "-e",
    (
        "rmarkdown::render("
        "input = commandArgs(trailingOnly=TRUE)[1], "
        "params = list("
        "patterns = commandArgs(trailingOnly=TRUE)[2], "
        "metabolome = commandArgs(trailingOnly=TRUE)[3]"
        "), "
        "output_file = commandArgs(trailingOnly=TRUE)[4])"
    ),
    rmd_path,
    patterns_path,
    metabolome_path,
    output_file_path,
]

# Execute the R command
subprocess.run(r_command, check=True)