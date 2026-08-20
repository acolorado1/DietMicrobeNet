import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Generate R Markdown report with parameters.")
parser.add_argument("--food_file", default="default.csv", help="Path to the data file.")
parser.add_argument("--output", default="filepath.html", help="Path to HTML report file")
args = parser.parse_args()

# Locate the Rmd file relative to this script's own location.
this_dir = os.path.dirname(os.path.abspath(__file__))
rmd_path = os.path.join(this_dir, "CompoundAnalysis_FooDB.Rmd")

# IMPORTANT: rmarkdown::render() knits with the working directory set to
# the .Rmd file's own directory (this_dir), NOT the directory the user ran
# this script from. If the .Rmd reads params$food_file with read_csv() (or
# similar) using relative-path semantics, we must resolve the user's input
# to an absolute path *before* passing it to R - anchored to the user's
# current working directory (where they actually ran this script), not to
# this_dir. Otherwise R will look for the file in the wrong place.
food_file_path = os.path.abspath(args.food_file)
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
        "params = list(food_file = commandArgs(trailingOnly=TRUE)[2]), "
        "output_file = commandArgs(trailingOnly=TRUE)[3])"
    ),
    rmd_path,
    food_file_path,
    output_path,
]

# Execute the R command
subprocess.run(r_command, check=True)