# Psychometric
These are the scripts for processing data for the psychometric grant

The scripts should be run in the following order
1. convertDicom
2. JC_preprocessing
3. JC_normc123
4. JC_FCpreproSetup -> JC_FCprepro
5. JC_FCCalcSetup -> JC_FCcalc
6. JC_ROIFCcalc
7. JC_3DMC.m
7. JC_rext_plot.m

Script #4 will call the JC_FCprepro function
Script #5 will call the JC_FCcalc function
Script #7 is used for quality control. This script goes through each subject's motion correction files and identifies outliers and writes them to a file for inspection.
Script #8 does...        


A new script will be added in the near future that will provide information about data quality and motion correction estimates.

##    convertDicom (script 1)
This is a shell script that calls dcm2niix to convert the raw dicom data to
from the scanner to nii files, which are used by all of the scripts below.


##    JC_preprocessing (script 2)
This script uses SPM12 to perform the following steps of preprocessing on all subjects

* if functional data doesn't have 42 slices, the subject will be skipped
* functional data is converted from integer to float32 datatype (to prevent precision less)
* the structural file is smoothed (with a 12mm FWHM Gaussian kernel)
* the smoothed version is then coregistered to the original (old) SPM12 T1 template (pulling along the structural scan)
* the smoothed version is then deleted (no longer needed)
* the structural file is then segmented (spm.spatial.preproc) using 6 tissue classes
* and the structural file is then normalized to MNI/ICBM space with a 1mm isotropic resolution
* the structural file (in original space) is then skull-stripped using imcalc and the c1/c2/c3 images as mask (where c1+c2+c3 >= 0.5)
* the functional data is then slice-timing corrected (order 2:2:42, 1:2:42; interleaved bottom-up, even-first)
* then the functional data is realigned (motion corrected) without reslicing, realign-to-mean, +write mean image
* the mean functional image is smoothed (also with a 12mm FWHM Gaussian kernel)
* and this smoothed version is then coregistered to the EPI template, pulling along all functional images
* next, the functional mean image is coregistered to the skull-stripped T1 file (still in subject space!), pulling along the functional data
* then the functional data is normalized to MNI/ICBM space with a 3mm isotropic resolution
* and finally, the normalized functional data is smoothed using a 6mm FWHM Gaussian kernel

##    JC_normc123 (script 3)
 Jochen wrote a second script to only normalize the c1/c2/c3 files to then perform the FCprepro step (which uses the GM/WM/CSF files *in the functional space* as global signal sources).

##    JC_FCpreproSetup (script 4)
This script finds the files used by the function: JC_FCprepro
### -->    JC_FCprepro (Uses):
JC_FCPREPRO  Preprocess functional data file(s) for functional connectivity
JC_FCPREPRO(PPFILE, RPFILE) removes the variance associated with
columns in RPFILE (text file or matfile) from PPFILE.

JC_FCPREPRO(PPFILE, RPFILE, TRFILT) also removes low-frequency drifts
present in the data up to TRFILT wavelength (using DCT filtering).

JC_FCPREPRO(PPFILE, RPFILE, TRFILT, GSFILES) also removes the mean
covariate from each file in GSFILES (resampled to the space in PPFILE)

This function uses SPM functions, and thus is dependent on SPM8 or
higher. Functions used include spm_vol, spm_read_vols, spm_write_vol,
spm_filter, and others.

The function will make a copy of (each of) the input file(s) with a
prefix of "fc" (e.g. "fcswraRUN.nii").

##    JC_FCCalcSetup (script 5)
This script finds the files used by the function JC_FCcalc
### -->JC_FCcalc (uses)
JC_FCPREPRO  Preprocess functional data file(s) for functional connectivity
JC_FCPREPRO(PPFILE, RPFILE) removes the variance associated with
columns in RPFILE (text file or matfile) from PPFILE.

JC_FCPREPRO(PPFILE, RPFILE, TRFILT) also removes low-frequency drifts
present in the data up to TRFILT wavelength (using DCT filtering).

JC_FCPREPRO(PPFILE, RPFILE, TRFILT, GSFILES) also removes the mean
covariate from each file in GSFILES (resampled to the space in PPFILE)

This function uses SPM functions, and thus is dependent on SPM8 or
higher. Functions used include spm_vol, spm_read_vols, spm_write_vol,
spm_filter, and others.

The function will make a copy of (each of) the input file(s) with a
prefix of "fc" (e.g. "fcswraRUN.nii").

6. JC_ROIFCcalc
7. JC_3DMC.m
7. JC_rext_plot.m

##    JC_3DMC.m (script 7)
USAGE:  This script gathers information about the motion correction estimates and then writes out the results to files.
