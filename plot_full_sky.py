import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import re
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from tqdm import tqdm
import scienceplots
plt.style.use(['science'])

# GBD
latitude = 13.603
longitude = 77.428
altitude = 694
location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg, height=altitude * u.m)

cache_dir = Path("cached_data")
cache_dir.mkdir(exist_ok=True)
cache_max_amps = cache_dir / "cache_max_amps.npy"
cache_ra = cache_dir / "cache_ra.npy"
cache_dec = cache_dir / "cache_dec.npy"

ra_bins = np.arange(0, 360 + 1.0, 1.0)
dec_bins = []
if cache_max_amps.exists() and cache_ra.exists() and cache_dec.exists():
    print("Cache files found. ")
    all_max_amps = np.load(cache_max_amps)
    all_ra = np.load(cache_ra)
    all_dec = np.load(cache_dec)
    print(f"Loaded {len(all_max_amps)} cached data points")

    elevation_dirs = list(Path(".").glob("elevation_*"))
    for elev_dir in elevation_dirs:
        elevation = float(elev_dir.name.split("_")[1])
        declination = 90.0 - elevation + latitude
        dec_bins.append(declination - 5.0)
        dec_bins.append(declination + 5.0)

else:

    print("Cache files not found. ")
    all_max_amps = []
    all_ra = []
    all_dec = []

    elevation_dirs = list(Path(".").glob("elevation_*"))
    print(f"Found {len(elevation_dirs)} elevation directories: {[d.name for d in elevation_dirs]}")
    for elev_dir in elevation_dirs:
        elevation = float(elev_dir.name.split("_")[1])
        declination = 90.0 - elevation + latitude
        dec_bins.append(declination - 5.0)
        dec_bins.append(declination + 5.0)
        print(f"Elevation: {elevation} degrees, Declination: {declination} degrees")
        trial_files = list(elev_dir.glob("trial_*.csv"))

        for trial_file in tqdm(trial_files, desc=f"Processing elevation {elevation}"):
            freqs = []
            col2 = []
            col3 = []
            amplitudes = []

            with open(trial_file, "r") as f:
                lines = f.readlines()

            data_start = False
            for line in lines:
                line = line.strip()

                if line.startswith("Frequency (MHz)"):
                    data_start = True
                    continue

                if not data_start or not line:
                    continue

                parts = line.split(",")
                if len(parts) < 3:
                    continue

                try:
                    fval = float(parts[0])
                    v2 = float(parts[1])
                    v3 = float(parts[2])
                except ValueError:
                    continue

                freqs.append(fval)
                col2.append(v2)
                col3.append(v3)
            amplitudes = [a - b for a, b in zip(col3, col2)]
            freqs = np.array(freqs)
            amplitudes = np.array(amplitudes)
            _, date_str, time_str = trial_file.stem.split("_")

            timestamp = f"{date_str} {time_str.replace('-', ':')}"
            obs_time = Time(timestamp)
            obs_time_utc = obs_time - 5.5 * u.hour

            if elevation <= 90:
                alt_val, az_val = elevation, 0
            else:
                alt_val, az_val = 180 - elevation, 180
            altaz = AltAz(alt=alt_val * u.deg, az=az_val * u.deg, location=location, obstime=obs_time_utc)

            skycoord = SkyCoord(altaz)
            ra = skycoord.icrs.ra.deg
            dec = skycoord.icrs.dec.deg

            # print(f"Trial: {trial_file.name}, Elevation: {elevation}, Dec: {dec:.2f}, RA: {ra:.2f}, Time (UTC): {obs_time_utc.iso}")

            maximum_amplitude = np.max(amplitudes)
            
            all_max_amps.append(maximum_amplitude)
            all_ra.append(ra)
            all_dec.append(dec)

    # Convert to numpy arrays and save to cache
    all_max_amps = np.array(all_max_amps)
    all_ra = np.array(all_ra)
    all_dec = np.array(all_dec)
    
    print("Saving computed data to cache...")
    np.save(cache_max_amps, all_max_amps)
    np.save(cache_ra, all_ra)
    np.save(cache_dec, all_dec)
    print(f"Cached {len(all_max_amps)} data points")

all_max_amps = np.array(all_max_amps)
all_ra = np.array(all_ra)
all_dec = np.array(all_dec)



dec_bins = np.unique(np.round(np.array(dec_bins), 1))
dec_bins.sort()

amp_sum, _, _ = np.histogram2d(
    all_dec,
    all_ra,
    bins=[dec_bins, ra_bins],
    weights=all_max_amps
)
counts, _, _ = np.histogram2d(
    all_dec,
    all_ra,
    bins=[dec_bins, ra_bins]
)
with np.errstate(invalid="ignore", divide="ignore"):
    amp_mean = amp_sum / counts

amp_mean[counts == 0] = np.nan

plt.figure(figsize=(12, 8))

mesh = plt.pcolormesh(ra_bins,dec_bins,amp_mean,cmap="turbo",shading="flat")

plt.colorbar(mesh, label="Peak Amplitude", orientation="horizontal", pad=0.1, fraction=0.05)
plt.xlabel("Right Ascension (degrees)", fontsize=20)
plt.ylabel("Declination (degrees)", fontsize=20)
plt.title("Sky Map", fontsize=24)

plt.gca().set_aspect("equal", adjustable="box")
plt.tight_layout()
plt.savefig("full_sky_map.png", dpi=300)
plt.show(block=False)
plt.pause(5)
plt.close()






