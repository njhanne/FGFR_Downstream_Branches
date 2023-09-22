# Golgi orientation
This code finds nuclei and their corresponding Golgi body, pairs them together, and calculates angle between their 
centroids. It's a bit more complicated than the proliferation code, but uses a lot of similar methods.

# Getting images and segmenting them with Cellpose
## Here are a number of steps I performed before this analysis
1. The images were captured using a Leica Stellaris confocal microscope. I have already converted the .lif files to tiffs 
and rotated them so that the bottom of the image is aligned to the FNP. Finally, I saved each channel separately using 
an ImageJ macro 'stack_channel_split.ijm' to speed up reading in the images to Cellpose. They are big files!
2. I compiled sample and image info into a .csv (GM130_sample_info.csv) that is used to tell the scripts what to do with each image. 
I think most of the columns are self-explanatory except:
   - z_space: this is the distance between slices in the z-direction taken by the confocal microscope.
   - max_distance: this tells the matching algorithm the cutoff distance between a nucleus and its Golgi body.
   - smallest_nuc: this tells the matching algorithm the cutoff for the smallest volume of a nucleus to consider.
   - angle_adjustment: this is an angle that will be subtracted from each Golgi angle measurement. The angle is the angle 
   of a line drawn between the left and right globular processes on a low resolution 'overview' image of the tissue section.
   See the example image below:
   ![Example of angle adjustment drawn on an overview image](/Readme_images/Golgi_baseline_angle_example.png)
3. Cellpose models are included for detecting nuclei and GM130. Both models were trained on images of Hoechst and GM130 
stained tissue sections from the chicken FNP. Every third image in the z-direction was labelled in ImageJ by hand. These 
labelled images were used to train with the Cellpose 'nuclei' model as the base. The models were refined further by 
'human-in-the-loop' training and correcting on unlabelled images. Using this method we only had to hand-correct ~5-10 images
before it was quite good at segmenting.
4. I will include the cellpose segmented images as well. I tried running the code without the GPU and it is very slow.
I imagine many people running this will not have access to a GPU.

## Actually starting the analysis
### 1. (optional) run_cellpose_Golgi.py
- Important! Don't forget to put 'False' instead of 'True' inside the code if you didn't install the cuda version of pytorch.
It's line 29 and says '### IMPORTANT! ###' above it. The images are quite large so the script will run for a pretty long time. 
Without CUDA it's going to take a very long time.
- This script will segment the Hoechst and GM130 labelled images and output mask files for each in the 'data/GM130/cellpose' directory.
The images have the same format as in the proliferation code – the pixel values of the image correspond to id of the segmented object.
So, in a nuclei image, pixel values of 0 are background, while non-zero values are different nuclei.

### 2. cellpose_angle_helpers.py
- There's alot happening in this script:
  1. First it loads in all the cellpose output mask files.
  2. Next it finds the centroids of each mask id - this is actually a pretty slow process, might even be slower than running
  cellpose! I've sped it up a bit by drawing bounding boxes around the mask so the entire image array isn't fed into the centroid
  finding code each time, but even finding those bounding boxes takes time.
  3. Match nuclei and Golgi centroids to their nearest neighbor. Sometimes the closest neighbor actually isn't the 'best' match.
  I get around this by using a [linear sum assignment](https://en.wikipedia.org/wiki/Assignment_problem) to find matching pairs, with
  the Euclidian distance as the 'cost' that the algorithm will minimize.
  4. Find a line between the paired nuclei and Golgi, convert the line to a unit vector, and use the i,j,k components to calculate the angle
  between them. It outputs all this data into a .csv for further analysis in R.

- Running this code should generate an image like the middle image below showing nuclei in blue and Golgi in yellow/green 
with a line drawn between the centroids.
![Diagram of centroids and filtering](/Readme_images/Golgi_filtering_example.png)

### 3. angle_analysis.R
- Again, I generally run R scripts line by line, but I tried to put most of the steps in functions to reduce clutter.
- First the code pulls in all the results .csv files generated from cellpose_angle_helpers.py, combines them into a large dataframe,
and cleans them up.
- I drew two rois in ImageJ to help filter the data:
  - The first is a freehand closed drawing around the mesenchyme, excluding
  any ectoderm or neurectoderm. In the right image in the above triptych, the red dots are cells excluded from this roi.
  - The second roi is a line I draw along the globular process. The R code will only analyze cells within 200um of this line, 
  which is represented as the black dots in the above image.
- Watson U-squared test is used to compare the distribution between the treated and contralateral sides of the face.
This analysis is only performed on the 2D data.
- The distribution of angles in 2D are plotted with a windrose histogram combining the relative number of pairs per sample
across treatment and side of the FNP.
- The distribution of angles in 3D are plotted with a Mollweide projection (a 'flatearth' plot) heatmap.

![windrose and flatearth plots](/Readme_images/Golgi_R_output.png)