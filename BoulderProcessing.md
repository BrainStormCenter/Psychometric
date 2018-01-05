Colorado notes and processing steps

# Processing steps
1. reorient/visualize (create a .mat file)
2. Art repair
3. INRIalign (creates: .mat, rp [movement], mean)
4. Normalize to epi template, spatial manipulation.
  * Use the SPM old normalize routine
  * the smoothing kernel used in this step influences edge detection etc..
  * The data written out during this step is NOT smoothed
5. smooth the normalized data



# Smoothing
  * Every time the data are resliced, the new written data are effectively smoothed.
  * The recommendation is to use a kernel that is 3-4 times the voxel size.


# Notes
  * when selecting files in spm you can use ^'prefix' to filter files
  * always be sure to set the select all the images in the .nii files. To do this set the range to '1:'number of volumes' (e.g., 1:300)
