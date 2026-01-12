"""


When naming a file of a planetary transit, make sure to name the object the name of the host star, otherwise you won't get the magnitude
because SIMBAD will check the planet, not the host star. Format must be WASP-XX, not WASP_XX.


"""

import os
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from datetime import datetime, timezone
from astropy import units as u

mjo_location = EarthLocation(lat=-43.9856 * u.deg, lon=170.4656 * u.deg, height=1029 * u.m) #I got the values from Wikipedia (Check)

FILTER_TO_FLUX = {'U': 'FLUX_U', 'B': 'FLUX_B', 'V': 'FLUX_V',
    'R': 'FLUX_R', 'I': 'FLUX_I',
    'J': 'FLUX_J', 'H': 'FLUX_H', 'K': 'FLUX_K',
    'g': 'FLUX_g', 'r': 'FLUX_r', 'i': 'FLUX_i', 'z': 'FLUX_z',}

def fits_scraper(fits_path):
    with fits.open(fits_path) as hdul:
        hdr = hdul[0].header
        object_name = hdr.get('OBJECT', None)
        filter_used = hdr.get('FILTER', None)
        if object_name is None:
            raise KeyError("OBJECT keyword not found in FITS header")
        if 'DATE-OBS' in hdr:
            obs_time = Time(hdr['DATE-OBS'], scale='utc')
        elif 'JD' in hdr:
            obs_time = Time(hdr['JD'], format='jd', scale='utc')
        else:
            raise KeyError("No DATE-OBS or JD found in FITS header")
    return object_name, obs_time, filter_used

def scrape_simbad(object_name, obs_time, filter_used):
    if obs_time is None:
        time_obs = Time(datetime.now(timezone.utc))
    elif isinstance(obs_time, Time):
        time_obs = obs_time
    else:
        time_obs = Time(obs_time, scale="utc")

    custom_simbad = Simbad()
    custom_simbad.TIMEOUT = 10
    custom_simbad.remove_votable_fields("coordinates")
    custom_simbad.add_votable_fields('ra(d)', 'dec(d)', 'otypes',
        'flux(U)', 'flux(B)', 'flux(V)', 'flux(R)', 'flux(I)',
        'flux(J)', 'flux(H)', 'flux(K)', 'flux(g)', 'flux(r)', 'flux(i)', 'flux(z)')

    try:
        result = custom_simbad.query_object(object_name)
        if result is None:
            raise ValueError(f"No SIMBAD entry found for {object_name}")

        ra_deg = result['RA_d'][0]
        dec_deg = result['DEC_d'][0]
        obj_type = result['OTYPES'][0].decode() if isinstance(result['OTYPES'][0], bytes) else result['OTYPES'][0]

        mag = None
        mag_band = '--'
        flux_field = FILTER_TO_FLUX.get(filter_used.lower())
        if flux_field in result.colnames:
            val = result[flux_field][0]
            if val is not None and not (hasattr(val, 'mask') and val.mask):
                mag = float(val)
                mag_band = filter_used.upper()

        if mag is None:
            for band in ['FLUX_V', 'FLUX_R', 'FLUX_I', 'FLUX_B', 'FLUX_U', 'FLUX_J', 'FLUX_H', 'FLUX_K']:
                if band in result.colnames:
                    val = result[band][0]
                    if val is not None and not (hasattr(val, 'mask') and val.mask):
                        mag = float(val)
                        mag_band = band[-1]
                        break

        if mag is None:
            mag = '--'
            mag_band = '--'

        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        altaz_frame = AltAz(obstime=time_obs, location=mjo_location) #Gets the azimuth and altitude
        altaz = coord.transform_to(altaz_frame)
        altitude = altaz.alt.deg
        azimuth = altaz.az.deg
        lst = time_obs.sidereal_time('apparent', longitude=mjo_location.lon)
        ha = (lst - coord.ra).wrap_at(24 * u.hour)
        airmass = altaz.secz


        diction = {
            'RA': coord.ra.to_string(unit=u.hour, sep=':'),
            'DEC': coord.dec.to_string(unit=u.deg, sep=':'),
            'EPOCH': 2000.0, #Hard coded in 
            'SIDTRACK': '--', #Hard code this in if needed
            'ALTITUDE': altitude,
            'AZIMUTH': azimuth,
            'LST': lst.to_string(unit=u.hour, sep=':'),
            'HA': ha.to_string(unit=u.hour, sep=':'),
            'AIRMSS': float(airmass.value) if hasattr(airmass, 'value') else airmass, #For some reason, the B&C header is AIRMSS, not AIRMASS
            'FLAGS': obj_type,
            'MAG': mag,
            'MAG_BAND': mag_band,
            # 'TELLIMIT': ['Limits'], 
            'UTTIME': time_obs.utc.iso
        }

    except Exception as e:
        print(f"SIMBAD error: {e}")
        diction = {'RA':'--','DEC':'--','EPOCH':'--','SIDTRACK':'--',
            'ALTITUDE':'--','AZIMUTH':'--','LST':'--','HA':'--',
            'AIRMSS':'--','FLAGS':'','MAG': '--', 'MAG_BAND':'--', 'UTTIME':'--'}
    return diction



def update_fits_header(fits_path, info_dict):
    with fits.open(fits_path, mode='update') as hdul:
        hdr = hdul[0].header
        for key, value in info_dict.items():
            if value is None:
                value = '--'
            hdr.set(key, value, "Value obtained from SIMBAD")
        hdul.flush()  # saves changes



def file_path(path):
    fits_files = [os.path.join(path, file) for file in os.listdir(path) if (file.lower().endswith(('.fit', '.fits')) and not file.lower().startswith('master'))]
    print(f"Found {len(fits_files)} FITS files")
    for fits_file in fits_files:
        print(f"\nProcessing: {fits_file}")
        try:
            object_name, obs_time, filter_used = fits_scraper(fits_file)
            header_info = scrape_simbad(object_name, obs_time, filter_used)
            update_fits_header(fits_file, header_info)
            print(f"{fits_file} updated successfully")
        except Exception as e:
            print(f"{fits_file} skipped because of {e}")



input_path = r"C:\Users\lmf53\Mount_John"
file_path(input_path)


# fits_file = r"C:\Users\lmf53\OneDrive - University of Canterbury\2026\Twin_Eyes\Reduction\20260106\BC_CMa-0002_i_30.fit"
# object_name, obs_time, filter_used = fits_scraper(fits_file)
# header_info = scrape_simbad(object_name=object_name, obs_time=obs_time, filter_used=filter_used)
# update_fits_header(fits_file, header_info)

# # for categories, values in header_info.items():
# #     print(f"{categories}={values}")

