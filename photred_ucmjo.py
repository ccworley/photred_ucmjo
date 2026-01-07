from os import killpg

import yaml
import copy
import matplotlib.pyplot as plt # A plotting library
import numpy as np # Numerical Python, great for vectorised equations
import warnings # To ignore our problems
from astropy.io import fits # Astronomy Python for opening ".fits" files
from astropy.time import Time
from photutils.background import Background2D, MedianBackground # Fitting background surfaces to astronomical images
from pathlib import Path
from collections import defaultdict
from collections.abc import Mapping


warnings.filterwarnings('ignore', category=RuntimeWarning) #Ignores some warnings

#%matplotlib widget

# ------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------
def load_calib_yaml(localdir):
    calib_yamlfile = Path(localdir + "/calibration_index.yaml")
    with open(calib_yamlfile, "r") as f:
        calib_yaml = yaml.safe_load(f)

    calib_yaml = copy.deepcopy(calib_yaml)  # protects against overwriting

    # print(calib_yaml.keys())
    # print(calib_yaml.get("flat", {}).keys())

    return calib_yaml


def generate_master_darks(calib_yaml, localdir, obsdate):

    dict_md = {}
    for imagetype, filt_dict in calib_yaml.items():
        print(imagetype)
        if 'dark' in imagetype:
            for filt, exp_dict in filt_dict.items():
                print(filt)
                for exptime, calib_entries in exp_dict.items():
                    print(exptime, calib_entries)

                    darks_data = []
                    for calib in calib_entries:
                        print(calib['file'])

                        hdul = fits.open(calib['file'])
                        hdr = hdul[0].header
                        darks_data.append(hdul[0].data.copy())
                        hdul.close()

                    # Take median of the darks data to create the master dark for the flats
                    master_dark = np.median(darks_data, axis=0)  # Median combine the darks

                    # Save master dark
                    md_filename = 'masterdark_' + str(exptime) + '.fits'
                    md_hdu = fits.PrimaryHDU(data=master_dark, header=hdr)
                    md_hdulist = fits.HDUList([md_hdu])  # List format (only 1 hdu in list)
                    md_hdulist[0].header['COMMENT'] = 'Master Dark: ' + str(exptime)
                    md_hdulist.writeto(localdir + "/" + obsdate + '/' + md_filename, overwrite=True)
                    md_hdulist.close()

                    # Save filename connect with exposure time for calling later for master flats and science
                    dict_md[exptime] = localdir + "/" + obsdate + '/' + md_filename

                    plt.figure()
                    plt.imshow(master_dark, cmap='gray', vmin=np.nanpercentile(master_dark, 2),
                               vmax=np.nanpercentile(master_dark, 98))  # 2d plotting
                    # plt.axis('off') # turning off axis labels
                    plt.title('Master Dark: ' + str(exptime))
                    plt.colorbar()
                    plt.savefig(localdir + "/" + obsdate + '/masterdark_' + str(exptime) + '.png')
                    #plt.show()
                    plt.close()

    return dict_md

def generate_master_flats(calib_yaml, localdir, obsdate, dict_md):

    dict_mf = {}
    for imagetype, filt_dict in calib_yaml.items():
        print(imagetype)
        if 'flat' in imagetype:
            for filt, exp_dict in filt_dict.items():
                print(filt)
                for exptime, calib_entries in exp_dict.items():
                    print(exptime, calib_entries)

                    #Find master dark that matches in exposure time
                    try:
                        md_file = dict_md[exptime]
                        hdrstr = 'ok'
                    except KeyError:
                        # Find closest exptime in available darks
                        md_exp = np.array(list(dict_md.keys()))
                        print(md_exp)
                        ediff = np.abs(exptime-md_exp)
                        mn_ediff = np.min(np.abs(exptime-md_exp))
                        alt_exp = [md_exp[l] for l,m in enumerate(ediff) if m == mn_ediff]
                        print(mn_ediff)
                        md_file = dict_md[alt_exp[0]]
                        hdrstr = 'Master Dark used had exptime=' + str(alt_exp)


                    md_hdul = fits.open(md_file)
                    md_hdr = md_hdul[0].header
                    md_data = md_hdul[0].data.copy()
                    md_hdul.close()

                    flats_data = []
                    for calib in calib_entries:
                        print(calib['file'])

                        hdul = fits.open(calib['file'])
                        hdr = hdul[0].header
                        flat_data = hdul[0].data.copy()
                        hdul.close()

                        # Ensure dataype is float  and subtract master dark with matching exposure time
                        flat_data = flat_data.astype(float) - md_data

                        print(np.nanmedian(flat_data))

                        # Check that the median count is optimal: between 24000 and 45000
                        if (np.nanmedian(flat_data) > 15000) & (np.nanmedian(flat_data) < 45000):
                            # print(np.nanmedian(flat_MasterDark))
                            flat_data /= np.nanmedian(flat_data)  # Normalise
                            flats_data.append(flat_data)


                    # Take median of the darks data to create the master dark for the flats
                    master_flat = np.median(flats_data, axis=0)  # Median combine the darks

                    # Save master flat
                    mf_filename = 'masterflat_' + filt + '_' + str(exptime) + '.fits'
                    mf_hdu = fits.PrimaryHDU(data=master_flat, header=hdr)
                    mf_hdulist = fits.HDUList([mf_hdu])  # List format (only 1 hdu in list)
                    mf_hdulist[0].header['COMMENT'] = 'Master Flat: '+ filt + '_' + str(exptime)
                    if 'ok' not in hdrstr:
                        mf_hdulist[0].header['COMMENT'] = hdrstr
                    mf_hdulist.writeto(localdir + "/" + obsdate + '/' + mf_filename, overwrite=True)
                    mf_hdulist.close()

                    # Save filename connect with exposure time for calling later for master flats and science
                    dict_mf[filt] = {'exptime':exptime, 'file':localdir + "/" + obsdate + '/' + mf_filename}

                    plt.figure()
                    plt.imshow(master_flat, cmap='gray', vmin=np.nanpercentile(master_flat, 2),
                               vmax=np.nanpercentile(master_flat, 98))  # 2d plotting
                    # plt.axis('off') # turning off axis labels
                    plt.title('Master Flat: '+ filt + '_' + str(exptime))
                    plt.colorbar()
                    plt.savefig(localdir + "/" + obsdate + '/masterflat_' + filt + '_' + str(exptime) + '.png')
                    #plt.show()
                    plt.close()

    return dict_mf


def get_science(all_files, HEADER_KEYS):
    global hdul, hdr, obj, filt, exptime
    for fits_path in all_files:
        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header

        # --- Identify calibration vs science ---
        imagetyp = hdr.get("IMAGETYP", "").lower()
        if "dark" in imagetyp or "flat" in imagetyp or "bias" in imagetyp or 'reduced' in str(fits_path):
            continue  # skip calibrations cleanly

        # --- Normalise OBJECT name ---
        obj_raw = hdr.get("OBJECT", "UNKNOWN")
        obj = obj_raw.strip().upper().replace(" ", "_")

        # --- Extract configuration ---
        cal_type = imagetyp.split()[0] if imagetyp else "science"
        filt = hdr.get(HEADER_KEYS["filter"], "NONE")
        exptime = hdr.get(HEADER_KEYS["exptime"], None)
        # print(obj, cal_type, filt, exptime)

        if exptime is None:
            continue  # science frames without EXPTIME are not useful

        # --- Time handling ---
        date_raw = hdr.get(HEADER_KEYS["date_obs"], None)
        if date_raw is not None:
            t = Time(date_raw, format="isot", scale="utc")
            date_obs = t.isot
            mjd = float(t.mjd)
        else:
            date_obs = None
            mjd = None

        # --- Entry ---
        entry = {
            "file": str(fits_path),
            "date_obs": date_obs,
            "mjd": mjd,
        }

        yaml_data[obj][cal_type][filt][float(exptime)].append(entry)

    # ------------------------------------------------------------------
    # Write YAML
    # ------------------------------------------------------------------
    with open(science_yamlfile, "w") as f:
        yaml.safe_dump(
            to_dict(yaml_data),
            f,
            sort_keys=True,
            default_flow_style=False,
        )

    return yaml_data


def generate_calibrated_science(science_yaml, localdir, obsdate, dict_md, dict_mf):

    dict_cs = {}
    for obj, imagetype_dict in science_yaml.items():

        for imagetype, filt_dict in imagetype_dict.items():
            print(imagetype)

            for filt, exp_dict in filt_dict.items():
                print(filt)

                # Find master flat that matches science in filter
                mf_hdrstr, mf_data = find_masterflat(dict_mf, filt)

                for exptime, science_entries in exp_dict.items():
                    print(exptime, science_entries)

                    #Find master dark that matches science in exposure time
                    md_hdrstr, md_data = find_masterdark(dict_md, exptime)

                    for science_entry in science_entries:
                        science_file = str(science_entry['file'])
                        print(science_file)

                        hdul = fits.open(science_file)
                        hdr = hdul[0].header
                        raw_data = hdul[0].data.copy().astype(float)

                        data = raw_data.copy() - md_data
                        data /= mf_data

                        bkg_estimator = Background2D(data, (50, 50), filter_size=(5, 5), bkg_estimator=MedianBackground())
                        median_background = bkg_estimator.background_median

                        data -= median_background

                        # # Save Calibrated Science Image
                        cs_filename = science_file.replace('.fit', '_reduced.fit')
                        cs_hdu = fits.PrimaryHDU(data=data, header=hdr)
                        cs_hdulist = fits.HDUList([cs_hdu])  # List format (only 1 hdu in list)
                        cs_hdulist[0].header['COMMENT'] = 'Reduced Science: ' + obj + ' '+ filt + '_' + str(exptime)
                        if 'ok' not in md_hdrstr:
                            cs_hdulist[0].header['COMMENT'] = md_hdrstr

                        cs_hdulist.writeto(cs_filename, overwrite=True)
                        cs_hdulist.close()

                        # Save filename connect with exposure time for calling later for master flats and science
                        dict_cs[science_file] = cs_filename

                        plt.figure()
                        plt.subplot(121)
                        plt.imshow(raw_data, cmap='gray', vmin=np.nanpercentile(raw_data, 2),
                                   vmax=np.nanpercentile(raw_data, 60))  # 2d plotting
                        # plt.axis('off') # turning off axis labels
                        plt.title('Raw Science: ' + obj + ' ' + filt + '_' + str(exptime),fontsize=7)
                        plt.colorbar(shrink=0.4)

                        plt.subplot(122)
                        plt.imshow(data, cmap='gray', vmin=np.nanpercentile(data, 2),
                                   vmax=np.nanpercentile(data, 60))  # 2d plotting
                        # plt.axis('off') # turning off axis labels
                        plt.title('Reduced Science: ' + obj + ' ' + filt + '_' + str(exptime),fontsize=7)
                        plt.colorbar(shrink=0.4)

                        plt.tight_layout()
                        plt.savefig(localdir + "/" + obsdate + '/reducedscience_'+ obj + '_' + filt + '_' + str(exptime) + '.png')
                        #plt.show()
                        plt.close()

    return dict_cs


def find_masterdark(dict_md, exptime) -> tuple[Any, str]:

    try:
        md_file = dict_md[exptime]
        hdrstr = 'ok'
    except KeyError:
        # Find closest exptime in available darks
        md_exp = np.array(list(dict_md.keys()))
        #print(md_exp)
        ediff = np.abs(exptime - md_exp)
        mn_ediff = np.min(np.abs(exptime - md_exp))
        alt_exp = [md_exp[l] for l, m in enumerate(ediff) if m == mn_ediff]
        #print(mn_ediff)
        md_file = dict_md[alt_exp[0]]
        hdrstr = 'Master Dark used had exptime=' + str(alt_exp)

    md_hdul = fits.open(md_file)
    #md_hdr = md_hdul[0].header
    md_data = md_hdul[0].data.copy()
    md_hdul.close()

    return hdrstr, md_data


def find_masterflat(dict_mf, filt) -> tuple[Any, str]:

    try:
        mf_file = dict_mf[filt]['file']
        hdrstr = 'ok'

    except KeyError:
        # If no flat in right filter can't reduce.
        print('Master Flat for filter ' + filt + ' not found')
        exit()

        # mf_exp = np.array(list(dict_mf.keys()))
        # # print(md_exp)
        # ediff = np.abs(exptime - mf_exp)
        # mn_ediff = np.min(np.abs(exptime - mf_exp))
        # alt_exp = [mf_exp[l] for l, m in enumerate(ediff) if m == mn_ediff]
        # # print(mn_ediff)
        # md_file = dict_md[alt_exp[0]]
        # hdrstr = 'Master Flat used had Filt=' + str(alt_exp)

    mf_hdul = fits.open(mf_file)
    # md_hdr = md_hdul[0].header
    mf_data = mf_hdul[0].data.copy()
    mf_hdul.close()

    return hdrstr, mf_data

# ------------------------------------------------------------------
# Nested dictionary factory: type -> filter -> exptime -> list
# ------------------------------------------------------------------

def exp_dict():
    return defaultdict(list)

def filter_dict():
    return defaultdict(exp_dict)

def type_dict():
    return defaultdict(filter_dict)

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

yaml_data = defaultdict(type_dict)

# ------------------------------------------------------------------
# Main
#
# Specify which OBSDATE for analysis
# If running locally, ensure files are located in similar structure to on tinsley
# Location of 'localdir' dictates all
# Use get_calibration_yaml.py to collate all possible flats and dark into a
# yaml file which this script interrogates
#
# Run script
# >python3 photred_ucmjo.py
# ------------------------------------------------------------------

obsdate = "20260105"
localdir = "/Users/ccw31/Documents/Data/UCMJO/Photometry/BC"  # BC -> octans
#localdir = "~/astro8/MJArchive/octans/"  # If on tinsley

science_yamlfile = Path(localdir + "/" + obsdate + "/science_index.yaml")

HEADER_KEYS = {
    "calibration_type": "IMAGETYP",
    "object": "OBJECT",
    "filter": "FILTER",
    "exptime": "EXPTIME",
    "date_obs": "DATE-OBS",
    #"julian_date": "JD",
}


# ------------------------------------------------------------------
# Load current yaml file with dictionary of dark and flat calibration files.
# ------------------------------------------------------------------
calib_yaml = load_calib_yaml(localdir)

# ------------------------------------------------------------------
# Generate updated Master Files for Darks+Exposure
# ------------------------------------------------------------------
dict_md = generate_master_darks(calib_yaml, localdir, obsdate)

# ------------------------------------------------------------------
# Generate updated Master Files for Flats+Filter+Exposure using Master Dark with matching Exposure
# ------------------------------------------------------------------
dict_mf = generate_master_flats(calib_yaml, localdir, obsdate, dict_md)

# ------------------------------------------------------------------
# ------------------------------------------------------------------
# Reduce Science Files using matching Master Dark and Master Flat
# ------------------------------------------------------------------
# Collect up Science files in OBSDATE directory
# ------------------------------------------------------------------
obsdir = [Path(localdir + "/" + obsdate)]
all_files = []  # Directories may contain science exposures. Filter out at next step
for d in obsdir:
    print(d)
    all_files.extend(d.glob("*.fit*"))

# ------------------------------------------------------------------
# Extract metadata and collect science files
# ------------------------------------------------------------------
science_yaml = get_science(all_files, HEADER_KEYS)

# ------------------------------------------------------------------
# Generate Reduced Science Images (no background subtraction yet)
# ------------------------------------------------------------------
dict_cs = generate_calibrated_science(science_yaml, localdir, obsdate, dict_md, dict_mf)
