import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Generate HTML report of microbe originating compounds.")
parser.add_argument("--node_file", default="default.csv", help="Path to the node data file.")
parser.add_argument("--edge_file", default="default.csv", help="Path to the edge data file.")
parser.add_argument("--output", default="filepath.html", help="Path to HTML report file")
args = parser.parse_args()

# Locate the Rmd file relative to this script's own location.
this_dir = os.path.dirname(os.path.abspath(__file__))
rmd_path = os.path.join(this_dir, "CompoundAnalysis_Microbes.Rmd")

# IMPORTANT: rmarkdown::render() knits with the working directory set to
# the .Rmd file's own directory (this_dir), NOT the directory the user ran
# this script from. If the .Rmd reads params$node_file / params$edge_file
# using relative-path semantics, we must resolve the user's inputs to
# absolute paths *before* passing them to R - anchored to the user's
# current working directory (where they actually ran this script), not to
# this_dir. Otherwise R will look for these files in the wrong place.
node_file_path = os.path.abspath(args.node_file)
edge_file_path = os.path.abspath(args.edge_file)
output_path = os.path.abspath(args.output)

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
        "node_file = commandArgs(trailingOnly=TRUE)[2], "
        "edge_file = commandArgs(trailingOnly=TRUE)[3]"
        "), "
        "output_file = commandArgs(trailingOnly=TRUE)[4])"
    ),
    rmd_path,
    node_file_path,
    edge_file_path,
    output_path,
]

# Execute the R command
subprocess.run(r_command, check=True)