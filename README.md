# Psychometric
Scripts for processing data for the psychometric grant
The scripts are to be run in the following order:
1. JC_preprocessing
2. JC_normc123
3. JC_FCprepro.m
4. JC_FCcalc.m
5. JC_ROIFCcalc.m
6. JC_rext_plot.m

A new script will be added in the near future that will provide information about data quality and motion correction estimates.

## JC_preprocessing
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
