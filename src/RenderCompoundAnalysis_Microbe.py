import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Generate HTML report of microbe originating compounds.")
parser.add_argument("--node_file", default="default.csv", help="Path to the node data file.")
parser.add_argument("--edge_file", default="default.csv", help="Path to the edge data file.")
parser.add_argument("--output", default="filepath.html", help="Path to HTML report file")
args = parser.parse_args()

this_dir = os.path.dirname(__file__)
rmd_path = os.path.join(this_dir, "CompoundAnalysis_Microbes.Rmd")
# Ensure output path is absolute
output_path = os.path.abspath(args.output)

# Construct the R command with parameters
r_command = f'''
Rscript -e "rmarkdown::render(
  input = '{rmd_path}',
  params = list(node_file = '{args.node_file}', edge_file = '{args.edge_file}'),
  output_file = '{output_path}')"
'''

# Execute the R command
subprocess.run(r_command, shell=True)
