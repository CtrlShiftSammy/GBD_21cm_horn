import re
from pathlib import Path

# input_file = "trials_results_2025-12-18 started at 10.08.39.csv"
input_files = list(Path(".").glob("trials_results_*"))
output_dir = Path("TRIALS")
output_dir.mkdir(exist_ok=True)

for input_file in input_files:
    with open(input_file, "r") as f:
        text = f.read()

    chunks = re.split(r"\n\s*\n\s*\n(?=Trial Number,\d+,)", text.strip())

    for i, chunk in enumerate(chunks, start=1):
        m = re.search(r"Trial Number,(\d+),", chunk)
        t = re.search(r"Observation Time,([^,]+),", chunk)
        obs_time = t.group(1).strip().replace(" ", "_").replace(":", "-")
        trial_id = m.group(1) if m else f"{i}"

        out_file = output_dir / f"trial_{obs_time}.csv"
        with open(out_file, "w") as f:
            f.write(chunk.strip() + "\n")

        print(f"Written {out_file}")
