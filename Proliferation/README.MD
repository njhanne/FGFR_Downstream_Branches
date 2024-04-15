# Cellular proliferation
This code finds cells in fluorescent IHC images and counts them. It's not all too complicated but there are some interesting
functions.

# Preparing images and segmenting them with Cellpose
## Here are a number of steps I performed before this analysis
1. The images were captured using the Leica fluorescent microscope in our lab. It is not confocal, the images were all captured
with the 40x objective. The fluorescent images are RGB since they are captured with a camera, not by the scope itself.
The pHH3 images are brightfield. The rgb_to_stack.py script converts all the RGB images to 8bit images then combines them into
a stack. It also inverts the pHH3 images so that the positive cells appear bright and the background appears dark.
![convert rgb to 8bit stack](/Readme_images/phh3_rbg_to_stack.png)
2. I used the 'in the loop' hand training function of Cellpose 2 to train models to segment the nuclei and pHH3 stained images.
3. Use the 'cellpoce_cmd.txt' command to find all the cells in the stacks. You need to run it twice, once for DAPI/Hoechst
labeled nuclei then again for the pHH3 labeled nuclei. This generates some pngs labelmaps where segmented object index
is the value of the pixel. So nuclei 221 will have a grey value of 221 in the image, the background will be 0.
![cellpose output](/Readme_images/proliferation_cellpose.png)

# Analyzing proliferation
4. The 'cellpose_proliferation_counting.py' creates a csv with the counts of nuclei and phh3 positive cells in each label 
image generated from Cellpose.
5. The results are analyzed in the 'proliferation.R' file. I think this is pretty self-explanatory.

# Other
There is another script in the Python directory called 'duplicate_check.py'. This compares all the RGB images to one another
to search for potential duplicates. It's rather brute force but not all that slow. I noticed two of the images seemed the
same and got paranoid... Turns out those two images were the only duplicates, but now we have a nifty script that can check for you.
It outputs a csv with all the comparisons and the confidence that the images are the same.