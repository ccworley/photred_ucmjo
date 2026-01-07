# Photometric Reduction for UCMJO Imaging Instruments

The photred_ucmjo.py is a Photometric reduction script for UCMJO.
The procedure has been adapted into a python script from the ASTR211 photometry jupyter notebooks.

It is designed to reduce raw images from the B&C and OC (not yet tested) into science images using flats and darks.

You can run it locally on your own computer or on tinsley (not yet tested).
if on your own computer, you need to set up the folders as if on tinsley.
The location of your 'localdir' should be equivalent to 'octans' in the MJArchive.

First, use the get_calibration_yaml.py to collate all possible flats and dark into a yaml file 
which this script interrogates. 

The calibration files can be spread over multiple OBSDATEs, or also in a Calibration_FLI/OBSDATE
folder as in octans. If pulling files from octans ensure files are put in matching OBSDATE locations.

The resulting calibration_index_yaml will be saved in the localdir. 
photred_ucmjo.py will look for it there.

photred_ucmjo.py currently runs on one OBSDATE but should do all objects and associated science files.
It will find the appropriate flats and darks as given in the calibration_index.yaml.

The output (masters and reduced images) are saved back into the OBSDATE folder with the science files.

To run script on command line
>python3 photred_ucmjo.py

>python3 get_calibration_yaml.py

The original development was using a .venv. Instructions to set up a .venv can/will be provided.

Development to do:
* Decide and update where output files will be saved and what they will be called
* Give option to search for most recent calibration files rather than any that match
* Further quality checks on the calibration files (only count limit on flats is currenlty implemented)
* Fill in missing header information (airmass, altitude,...)
* Apply WCS (see astr211 notebook in jupyternbs folder)
* Use a configuration file and/or command arguments for directories to look at rather than inserting in code
* ....other things?
