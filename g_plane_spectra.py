import math
import time
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, Galactic, get_body_barycentric_posvel, CartesianRepresentation
import astropy.units as u
from astropy.time import Time
from tqdm import tqdm
import scienceplots
plt.style.use(['science'])

# GBD location
latitude = 13.603
longitude = 77.428
altitude = 694
location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=altitude*u.m)
f0 = 1420.405
c = 299792.458  # km/s
def find_galactic_plane_crossings(dec):
    """Find the two RA values where galactic latitude = 0 at given declination"""
    ra_test = np.linspace(0, 360, 3600)
    coords = SkyCoord(ra=ra_test * u.deg, dec=dec * u.deg, frame='icrs')
    galactic_coords = coords.galactic
    glat = galactic_coords.b.deg
    
    # Find zero crossings
    sign_changes = np.where(np.diff(np.sign(glat)))[0]
    ra_crossings = []
    for idx in sign_changes:
        # Interpolate to find exact crossing
        ra1, ra2 = ra_test[idx], ra_test[idx + 1]
        glat1, glat2 = glat[idx], glat[idx + 1]
        ra_crossing = ra1 - glat1 * (ra2 - ra1) / (glat2 - glat1)
        ra_crossings.append(ra_crossing)
    
    return ra_crossings[:2] if len(ra_crossings) >= 2 else ra_crossings

def load_spectrum(trial_file):
    """Load frequency and amplitude from a trial file"""
    freqs = []
    col2 = []
    col3 = []
    
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
    return np.array(freqs), np.array(amplitudes)

def get_ra_from_filename(trial_file):
    """Extract RA from trial filename timestamp"""

    
    _, date_str, time_str = trial_file.stem.split("_")
    timestamp = f"{date_str} {time_str.replace('-', ':')}"
    obs_time = Time(timestamp)
    obs_time_utc = obs_time - 5.5 * u.hour
    
    elev_dir = trial_file.parent
    elevation = float(elev_dir.name.split("_")[1])
    
    location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg, height=altitude * u.m)
    
    if elevation <= 90:
        alt_val, az_val = elevation, 0
    else:
        alt_val, az_val = 180 - elevation, 180
    
    altaz = AltAz(alt=alt_val * u.deg, az=az_val * u.deg, location=location, obstime=obs_time_utc)
    skycoord = SkyCoord(altaz)
    ra = skycoord.icrs.ra.deg
    
    return ra

# Process each elevation directory
elevation_dirs = sorted(Path(".").glob("elevation_*"))
print(f"Found {len(elevation_dirs)} elevation directories")

# fig, axes = plt.subplots(len(elevation_dirs), 2, figsize=(14, 4 * len(elevation_dirs)))
# if len(elevation_dirs) == 1:
#     axes = axes.reshape(1, -1)

for idx, elev_dir in enumerate(elevation_dirs):
    
    elevation = float(elev_dir.name.split("_")[1])
    declination = 90.0 - elevation + latitude
    
    print(f"\nProcessing elevation {elevation}° (dec={declination:.2f}°)")
    
    # Find galactic plane crossings
    ra_crossings = find_galactic_plane_crossings(declination)
    print(f"Galactic plane crossings at RA: {ra_crossings[0]:.2f}°, {ra_crossings[1]:.2f}°")
    
    trial_files = list(elev_dir.glob("trial_*.csv"))
    
    # Process each crossing
    for crossing_idx, ra_crossing in enumerate(ra_crossings):
        spectra = []
        velocity = []
        velocity_lsr = []
        plt.figure(figsize=(12, 6))
        
        for trial_file in tqdm(trial_files, desc=f"Elev {elevation}°, crossing {crossing_idx + 1}"):
            ra = get_ra_from_filename(trial_file)
            
            # Check if within ±5° of crossing (handle wrap-around at 0/360)
            ra_diff = abs(ra - ra_crossing)
            if ra_diff > 180:
                ra_diff = 360 - ra_diff
            
            if ra_diff <= 2.50:
                freqs, amps = load_spectrum(trial_file)
                if len(freqs) > 0:
                    spectra.append(amps)
                    velocity.append((c * (f0 - freqs) / f0))

                    print(f"Trial file name: {trial_file}")
                    timestamp = f"{trial_file.stem.split('_')[1]} {trial_file.stem.split('_')[2].replace('-', ':')}"
                    obs_time = Time(timestamp)
                    print(f"Observation Time (UTC): {obs_time - 5.5 * u.hour}")
                    gcrs = location.get_gcrs(obs_time - 5.5 * u.hour)
                    gcrs_icrs = gcrs.transform_to(ICRS())
                    v_obs_geocentric = gcrs.velocity.to_cartesian()
                    earth_posvel = get_body_barycentric_posvel('earth', obs_time - 5.5 * u.hour)
                    v_earth_bary = earth_posvel[1].xyz.to(u.km/u.s)
                    v_total_bary = v_obs_geocentric.xyz.to(u.km/u.s) + v_earth_bary 
                    target = SkyCoord(ra=ra * u.deg, dec=declination * u.deg, frame='icrs')
                    target_ra = math.radians(target.icrs.ra.deg)
                    target_dec = math.radians(target.icrs.dec.deg)
                    los_unit = target.cartesian.xyz / target.cartesian.norm()
                    

                    v_bary_along_los = np.dot(v_total_bary.to(u.km/u.s).value, los_unit.value) * u.km/u.s
                    # sun's peculiar velocity wrt LSR (Schönrich et al. 2009)
                    v_sun_mag = 18.044 * u.km/u.s
                    sun_apex = SkyCoord(ra=267.058285*u.deg, dec=23.096545*u.deg, frame='icrs')
                    sun_unit = sun_apex.cartesian.xyz / sun_apex.cartesian.norm()

                    v_sun_along_los = v_sun_mag * np.dot(sun_unit.to_value(u.one), los_unit.to_value(u.one))
                    # print(f"v_sun_along_los: {v_sun_along_los}")
                    print(f"v_bary_along_los: {v_bary_along_los}")
                    print(f"v_sun_along_los: {v_sun_along_los}")
                    v_corr_lsr = v_bary_along_los + v_sun_along_los
                    print(f"Correction term: {v_corr_lsr:.2f}")
                    velocity_lsr.append(velocity + float(v_corr_lsr.to(u.km/u.s).value))

        if len(spectra) > 0:
            # all spectra
            spectra = np.array(spectra)
            velocity_lsr = np.array(velocity_lsr)
            for v_lsr, spec in zip(velocity_lsr, spectra):
                plt.plot(v_lsr, spec, color='gray', alpha=0.3, linewidth=0.5)
            # Average spectra
            avg_spectrum = np.mean(spectra, axis=0)
            coordinate = SkyCoord(ra=ra_crossing * u.deg, dec=declination * u.deg, frame='icrs')
            galactic_lat = coordinate.galactic.b.deg
            galactic_lon = coordinate.galactic.l.deg
            np.savez(f"galactic_plane_spectra_{ra_crossing:.1f}_{declination:.1f}.npz", 
                     freqs=freqs, avg_spectrum=avg_spectrum, spectra=spectra, 
                     ra=ra_crossing, declination=declination, glatitude=galactic_lat, glongitude=galactic_lon, elevation=elevation)
            
            # plt.plot(freqs, avg_spectrum, linewidth=2, color='navy')
            plt.xlabel("Frequency (MHz)", fontsize=24)
            plt.ylabel("Amplitude", fontsize=24)
            plt.title(rf"Elev {elevation}° | RA={ra_crossing:.1f}° | Dec={declination:.1f}° | $\ell$={galactic_lon:.1f}° | $b$={galactic_lat:.1f}° | (n={len(spectra)})", fontsize=24)
            plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"galactic_plane_spectra_{ra_crossing:.1f}_{declination:.1f}.png", dpi=300, bbox_inches='tight')
        print(f"\nSaved galactic_plane_spectra_{ra_crossing:.1f}_{declination:.1f}.png")
        plt.show(block=False)
        plt.pause(10)
        plt.close()
        time.sleep(0.5)
