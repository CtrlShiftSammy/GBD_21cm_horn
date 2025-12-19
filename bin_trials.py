from pathlib import Path
from datetime import datetime, timedelta
import shutil


log_file = Path("log.txt")
elevation_bins = []

with open(log_file, "r") as f:
    lines = f.readlines()[1:]  
    for line in lines:
        parts = line.strip().split(",")
        if len(parts) >= 4:
            date_str = parts[1]  
            time_str = parts[2]  
            elevation = parts[3]
            
            
            dt = datetime.strptime(f"{date_str} {time_str}", "%d-%m-%y %H:%M")
            end_time = dt + timedelta(hours=24)
            
            elevation_bins.append({
                "elevation": elevation,
                "start": dt,
                "end": end_time
            })


trials_dir = Path("TRIALS")
trial_files = list(trials_dir.glob("trial_*.csv"))

for trial_file in trial_files:
    filename = trial_file.stem  
    time_part = filename.replace("trial_", "")
    
    try:
        trial_time = datetime.strptime(time_part, "%Y-%m-%d_%H-%M-%S")
        matched_elevation = None
        for bin_info in elevation_bins:
            if bin_info["start"] <= trial_time < bin_info["end"]:
                matched_elevation = bin_info["elevation"]
                break
        
        if matched_elevation:
            
            elev_folder = Path(f"elevation_{matched_elevation}")
            elev_folder.mkdir(exist_ok=True)
            
            
            dest_file = elev_folder / trial_file.name
            shutil.move(str(trial_file), str(dest_file))
            print(f"Moved {trial_file.name} to elevation_{matched_elevation}/")
        else:
            print(f"No matching elevation found for {trial_file.name} (time: {trial_time})")
    
    except ValueError as e:
        print(f"Could not parse time from {trial_file.name}: {e}")

print("\nBinning complete!")
