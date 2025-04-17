# Cellular proliferation
This code finds cells in fluorescent IHC images and counts them. It's not complicated but there are some interesting
functions. It is almost the same as the proliferation code.

# Preparing images and segmenting them with Cellpose
## Here are a number of steps I performed before this analysis
1. The images were captured using the Leica fluorescent microscope in our lab. It is not confocal, the images were all captured
with the 40x objective. The fluorescent images are RGB since they are captured with a camera, not by the scope itself.
In this case we have two sets of images, one with Hoecsht and TUNEL, and one with TUNEL alone. The rgb to stack image subtracts
the TUNEL only image from the combined image, converts everything to 8bit, then converts to a stack.

The bottom left image is the combined RGB, the middle is the TUNEL only RGB, and the right is the 8bit stack with the TUNEL 
channel turned off.
![convert rgb to 8bit stack](/Readme_images/TUNEL_rgb_to_stack.png)
2. I used the 'in the loop' hand training function of Cellpose 3 to train models to segment the nuclei and pHH3 labeled images.
I did not train a model for the TUNEL, I just used the same pHH3 model which works fine for segmenting TUNEL, too.
3. Use the 'cellpose_cmd.txt' command to find all the cells in the stacks. You need to run it twice, once for DAPI/Hoechst
labeled nuclei then again for the TUNEL labeled nuclei. It is easier to just use the 'run_cellpose_proliferation.py' instead, but
note you will need to edit the script to specify if you are using cuda or not.

This generates some pngs labelmaps where segmented object index
is the value of the pixel. So nuclei 221 will have a grey value of 221 in the image, the background will be 0. The bottom 
left image is the labelmap for nuclei and the bottom right is for the TUNEL. The yellow rois in both images show the 
included outer area and excluded bead. Note that the high amount of apoptosis in the area where the globular process is 
undergoing fusion with the maxillary process has not been included in this analysis, since this apoptosis is likely not
due to the bead but normal developmental programmed cell death.
![cellpose output](/Readme_images/apoptosis_cellpose.png)

# Analyzing proliferation
4. As seen in the images above, I draw one or two rois in each image. The first outlines the area where apoptosis is analyzed.
The ectoderm and areas of fusion are excluded. The second roi was drawn around the bead if it was present in the image, which
is excluded from analysis. See the yellow rois in the above diptych.
5. The results are analyzed in the 'apoptosis.R' file, which is annotated.