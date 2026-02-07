import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description="Generate HTML report of microbe originating compounds.")
parser.add_argument("--patterns", default="default.csv", help="Path to the graph_results.csv")
parser.add_argument("--metabolome", default="default.csv", help="Path to CSV file containing one column of KEGG compounds.")
parser.add_argument("--output", default="filepath.html", help="Path to HTML report file")
args = parser.parse_args()

this_dir = os.path.dirname(__file__)
rmd_path = os.path.join(this_dir, "MetabolomeComparison_Report.Rmd")
# Ensure output path is absolute
output_path = os.path.abspath(args.output)

# Construct the R command with parameters
r_command = f'''
Rscript -e "rmarkdown::render(
  input = '{rmd_path}',
  params = list(patterns = '{args.patterns}', metabolome = '{args.metabolome}'),
  output_file = '{output_path}/MetabolomeComparison.html')"
'''

# Execute the R command
subprocess.run(r_command, shell=True)