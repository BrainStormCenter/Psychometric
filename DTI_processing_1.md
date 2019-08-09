# DTI Processing



		CREATED BY:	Huiling Peng
		CREATED ON:	2018-11-15

		USAGE:	To process DTI data
---

The following will provide information on how to get FA, AD, MD images from your DTI data. After you have these files ready, then you can decide what kind of data analysis you want and start from there!!

1. I use dcm2nii to convert all raw ima or dcm files to \*.nii.gz. If you choose in your protocol to output FA, MD, etc, you will have many of nii.gz files. The only one we need is the 4D raw DWI nii.gz file. It has the biggest file size. dcm2nii will also extract \*.bval and \*.bvec files. These three files are necessary for DTI processing

2. I would recommend to rename these three files.  
a. rename the 4D dwi nii.gz file to: rawdata.nii.gz  
 b. rename the \*.bval to: SubXXDTI.bvals (note the 's')  
 c. rename the \*.bvec to SUBXXXDTI.bvecs

For the following data processing, if you are not familiar with FSL command line.  I would recommend to use the FSL GUI. Just type ‘fsl &’ to run it if it is correctly installed. Here assume we put the 3 input files mentioned above under a directory named ‘DTI’. I will use DTI/ to represent the absolute path for the following description.

1. Eddy current correction.

1)Press the "FDT Diffusion" button from the main FSL menu to show the diffusion processing window
2) select EDDY current correction
3) select your ‘DTI/rawdata’ as input, and type ‘DTI/data’ as corrected output name.
4) The reference volume is typically 0
5) Press ‘Go” to run the eddy current correction and wait for it to be compete.

2. brain extraction. This ensures that we only compute tensors inside the brain.

1) Press ‘BET’ button from the main FSL window.
2) choose /DTI/data as input
3) type ‘/DTI/nodif_brain’ as output
4) Under Advanced options, check ‘Output binary brain mask image’
5) The Fractional intensity threshold should be set to 0.3 (we would rather including scalp than removing part of brain)
6) Press ‘Go’ and wait for it to complete. You will have DTI/nodif_brain and DTI/nodif_brain_mask created.

3. diffusion tensor reconstruction

1) make sure you have data.nii.gz, bvecs, bvals, nodif_brain_mask under your DTI/ directory
2) Press the “FDT Diffusion” button from the main FSL menu and select ‘DTIFIT Reconstruct diffusion tensor’.
3) Select the input directory ‘DTI/’
4 Press ‘Go’ and wait for it to be complete.

Check you DTI/ directory, you should have all your output of dti_FA.nii.gz dti_S0.nii.gz, dti_V1.nii.gz etc. You want to check your images, especially FA and V1 using fslview.

I copied from the website  http://www.cabiatl.com/Resources/Course/tutorial/html/dti.html how to show these two images.





 Vieiwing DTI data with FSLview
FSL has written the results as NIfTI format files, and you can open then in FSLView.
Your files should be saved in your data folder. It will save one file with the title 'FA', which refers to the Fractional Anistropy map. It also saves the three orthogonal vectors (V1, V2, V3).
From the Unix command line, change directory to the folder with your DTI data (e.g. 'cd ~/tutorial')
Launch FSL by typing 'FSL &' from the command line.
Press the FSLView from the main FSL menu to show launch FSLView.
Choose File/Open and select the FA map (e.g. DTI_FA).
Choose File/Add to add an overlay - select the V1 file (DTI_V1) to show principle vectors.
In box on the bottom, you can see what you have loaded (FA and V1); to turn the overlay off, just click on the eyeball. 
Select the info button (i with circle) after highlighting V1.
Here is how to show directions as lines
Image type: Diffusion Tensor
DTI Display options: LINES. 
Here is how to show the data as colors
Image type: Diffusion Tensor
DTI display options, Display as: RGB
DTI display options, Modulation: dti_FA  In this view, regions with little FA anistropy (e.g. CSF), whereas tissue with extreme FA are bright (fiber tracts). Color codes to give directions: Red is left to right. Green is front to back. Blue is head to foot.
