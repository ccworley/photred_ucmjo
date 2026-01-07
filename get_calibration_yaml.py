import yaml
from pathlib import Path
from astropy.io import fits
from astropy.time import Time
from collections import defaultdict
from collections.abc import Mapping

"""
Calibration Index Builder
=========================

This script scans one or more directories containing FITS calibration files
(e.g. bias, dark, flat frames) and builds a structured YAML index describing
their essential metadata.

The index is designed specifically to support *algorithmic selection* of
calibration files for science observations, rather than simple bookkeeping.

Key design principles
---------------------
1. Calibrations are grouped hierarchically by:
      calibration type  ->  filter  ->  exposure time

   This mirrors how calibrations are actually matched in practice and avoids
   expensive or error-prone scanning of irrelevant files.

2. Each calibration frame is indexed individually, not averaged or collapsed.
   Exact observation times are preserved so that calibration frames can be
   selected based on temporal proximity to a science exposure.

3. Observation times are stored both as ISO strings and as MJD values.
   MJD enables fast and numerically stable sorting by |Δt| when selecting the
   closest calibration frames in time.

4. Exposure time is used as a dictionary key rather than a free attribute.
   This allows deterministic matching (or tolerance-based matching) without
   ambiguity or string parsing.

Why this approach works well
----------------------------
• Enables reproducible, explainable calibration selection
• Scales naturally to large data volumes
• Separates *indexing* from *selection logic*
• Avoids premature creation of master calibrations
• Closely resembles professional CalDB-style systems

The resulting YAML file is intended to be consumed later by a calibration
selection routine that matches science frames on:
   - calibration type
   - filter
   - exposure time (exact or within tolerance)
   - temporal proximity (closest MJD)

This makes the calibration process explicit, auditable, and robust.
"""
# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------

localdir = "/Users/ccw31/Documents/Data/UCMJO/Photometry/BC"  # BC -> octans
#localdir = "~/astro8/MJArchive/octans/"  # If on tinsley

# OBSDATE and Calibration directories should be named the same and in the same structure as in MJArchive
calibration_dirs = [
    Path(localdir + "/20260105"),
    Path(localdir + "/Calibration_FLI/20251222"),
    Path(localdir + "/20251122")
    #Path("/data/calibration/flat"),
]

yaml_outfile = Path(localdir + "/calibration_index.yaml")

HEADER_KEYS = {
    "calibration_type": "IMAGETYP",
    "filter": "FILTER",
    "exptime": "EXPTIME",
    "date_obs": "DATE-OBS",
    #"julian_date": "JD",
}


# ------------------------------------------------------------------
# Nested dictionary factory: type -> filter -> exptime -> list
# ------------------------------------------------------------------

def exp_dict():
    return defaultdict(list)

def filter_dict():
    return defaultdict(exp_dict)

def to_dict(obj):
    """
    Recursively convert any mapping type (including defaultdict)
    into a plain dict so it can be safely serialized to YAML.
    """
    if isinstance(obj, Mapping):
        return {k: to_dict(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [to_dict(v) for v in obj]
    else:
        return obj

yaml_calib = defaultdict(filter_dict)


# ------------------------------------------------------------------
# Collect calibration files
# ------------------------------------------------------------------

all_files = []  # Directories may contain science exposures. Filter out at next step
for d in calibration_dirs:
    print(d)
    all_files.extend(d.glob("*.fit*"))

# ------------------------------------------------------------------
# Extract metadata
# ------------------------------------------------------------------

for fits_path in all_files:
    with fits.open(fits_path) as hdul:
        hdr = hdul[0].header

    cal_type = hdr.get("IMAGETYP", "UNKNOWN").split(' ')[0].lower()
    filt = hdr.get(HEADER_KEYS["filter"], "NONE")
    exptime = hdr.get(HEADER_KEYS["exptime"], None)
    #jd = hdr.get(HEADER_KEYS["julian_date"], None)
    print(cal_type, filt, exptime) #, jd)

    if 'dark' in cal_type or 'flat' in cal_type:
        date_raw = hdr.get(HEADER_KEYS["date_obs"], None)
        if date_raw is not None:
            t = Time(date_raw, format="isot", scale="utc")
            date_obs = t.isot
            mjd = float(t.mjd)
        else:
            date_obs = None
            mjd = None

        entry = {
            "file": str(fits_path),
            "date_obs": date_obs,
            "mjd": mjd,
        }

        yaml_calib[cal_type][filt][float(exptime)].append(entry)


# ------------------------------------------------------------------
# Write YAML
# ------------------------------------------------------------------

with open(yaml_outfile, "w") as f:
    yaml.safe_dump(
        to_dict(yaml_calib),
        f,
        sort_keys=True,
        default_flow_style=False,
    )

print(f"Wrote time-aware calibration index to {yaml_outfile}")
