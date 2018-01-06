Colorado notes and processing steps

# Processing steps
1. reorient/visualize (create a .mat file)
  * open the raw data with the 'display' button
  * put crosshair on AC and press 'Reorient'
2. Art repair
3. INRIalign (creates: .mat, rp [movement], mean)
  * Reslice only
  * sinc interpolation
  * create mean image only
4. Normalize to epi template, spatial manipulation.
  * Use the SPM old normalize routine
  * be sure to delete the original normalize module
  * Source image = mean images
  * images to write = functional data set
  * template = epi template
  * nonlinear frequency cutoff = 45 (recommended)
  * set voxel size to [3 3 3]
  * the smoothing kernel used in this step influences edge detection etc..
  * The data written out during this step is NOT smoothed
5. smooth the normalized data
6. Specify the first level analyses
  * set output directory
  * units for design = scans
  * interscan interval = TR
  * choose scans with the 'sw' prefix
  * onset = use 'spm_load' to select text file with stimulus onsets
  * see cheat_sheet_fMRImodelspecification.doc for additional details




# Smoothing
  * Every time the data are resliced, the new written data are effectively smoothed.
  * The recommendation is to use a kernel that is 3-4 times the voxel size.


# Choosing a basis set
  * There is a step that has to be done outside SPM, but I'm not sure what that is right now.
  * Tor recommends using a FIR model and Vince will talk about a 'boost' procedure.

# Notes
  * when selecting files in spm you can use ^'prefix' to filter files
  * always be sure to set the select all the images in the .nii files. To do this set the range to '1:'number of volumes' (e.g., 1:300)
  * Do NOT yolk stimuli to the TR
  * SPM will bin each volume into microtime units
  * Jittering stimulus onset within the TR can effectively increase sampling rate
  * high pass filtering remove low frequencies out of the data to increase the sensitivity to the BOLD signal.
  * ask about how to do the derivative outside of SPM.
  * Scaling is not needed and should really never been done.
  * in the contrast manager, make sure the contr
  * Presentation is the recommended stimulus presentation software



# Processing data
## auditory oddball paradigm
  * 12% target (1000 MHz 200ms tone)
  * 12% novel (random 700ms tone)
  * 76% standard (500ms tone)
