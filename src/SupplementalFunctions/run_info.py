import datetime
import json

def write_run_info(version, config, output):
    with open(output, "w") as f:
        f.write(f"version: {version}\n")
        f.write(f"date: {datetime.datetime.now()}\n")
        f.write("config:\n")
        f.write(json.dumps(config, indent=2))