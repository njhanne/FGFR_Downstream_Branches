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
stained tissue sections from the chick FNP. Every third image in the z-direction was labelled in ImageJ by hand. These 
labelled images were used to train with the Cellpose 'nuclei' model as the base. The models were refined further by 
'human-in-the-loop' training and correcting on unlabelled images. Using this method we only had to hand-correct ~5-10 images
before it was quite good at segmenting.
4. I will include the cellpose segmented images as well. I tried running the code without the GPU and it is very slow.
I imagine many people running this will not have access to a GPU, unfortunately.

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
  finding code each time, but even finding those bounding boxes takes time. Could probably just get center of bounding box
  to approximate centroid, but let's be accurate and get the true centroid.
  3. Match nuclei and Golgi centroids to their nearest neighbor. Sometimes the closest neighbor actually isn't the 'best' match.
  I get around this by using a [linear sum assignment](https://en.wikipedia.org/wiki/Assignment_problem) to find matching pairs, with
  the Euclidian distance as the 'cost' that the algorithm will minimize.
  4. Find a line between the paired nuclei and Golgi, convert the line to a unit vector, and use the i,j,k components to calculate the angle
  between them. It outputs all this data into a .csv for further analysis in R.

- Running this code will generate an image like the middle image below showing nuclei in blue and Golgi in yellow/green 
with a line drawn between the centroids.
![Diagram of centroids and filtering](/Readme_images/Golgi_filtering_example.png)

### 3. angle_analysis.R
This contains many functions and tools for analyzing the data, many of which are not going to appear in the manuscript,
but I do not have the heart to delete them. They may be useful to someone! My pipeline has changed a lot in the last few
months and the code has many deprecated sections, but again, I feel they could be useful to other scientists.
- I generally run R scripts line by line, but I tried to put most of the steps in functions to reduce clutter.
You can safely run the first ~1200 lines at once to load up all the libraries and functions.
- First the code pulls in all the results .csv files generated from cellpose_angle_helpers.py, combines them into a large dataframe,
and cleans them up.
- I drew two rois in ImageJ to help filter the data:
  - The first is a freehand closed drawing around the mesenchyme, excluding
  any ectoderm or neurectoderm. In the right image in the above triptych, the red dots are cells excluded from this roi.
  - DEPRECATED: The second roi is a line I draw along the globular process. The R code will only analyze cells within 200um of this line, 
  which is represented as the black dots in the above image.
- NEW (BETTER): I drew 4 pts in ImageJ on the low-res overview image of the entire FNP and in the same locations on the 
high-res confocal z-stacks. These points are used to transform the two z-stacks to a single coordinate system of the entire
tissue (see image below, you can see the outline of the two z-stacks in the overview image).
![overview image with landmarks and z-stack locations highlighted](/Readme_images/overview_angle_lms.png)
  - On the overview image I draw 8 line rois ('octiles') around the tissue edges. These rois are used to determine Golgi
  'positional angle', which will be described more later. These rois are used to filter nuclei a certain distance from 
  the FNP ectoderm, similar but less subjective than the previous roi method. Note in the below plot of nuclei positions,
  the black are 200 um from the baseline, and the pink? salmon? ones are further. The ectoderm and neurectoderm nuclei
  have already been removed before transformation so they do not appear in this graph (red dots in triptych above). A third 
  landmark would help anchor the 'y' scale better, but it is not necessary for this analysis.
![overview image with octile rois](/Readme_images/overview_octile_filtered.png)

- For most of the analysis the z-component of the angles are removed.
- Golgi angle is simply the angle formed between the nucleus and the Golgi. Positional orientation, or positional 'angle',
is easier to describe with an image than writing, so see the cartoon below. Basically it is the distance along the octile 
roi line converted into an angle. As such, it is not truly angular data, but just a convenient way to present the data.
The positional angle describes where along the tissue the Golgi is 'pointing'.
<img src="/Readme_images/positional_angle_diagram.png" width=50% height=50%>
- Watson U-squared test is used to compare Golgi angle distribution between the treated and contralateral sides of the face.
This analysis is only performed on the 2D data.
- The distribution of angles and positional angle in 2D are plotted with a windrose histogram combining the relative number 
of pairs per sample across treatment and side of the FNP.
- The distribution of Golgi angle in 3D are plotted with a Mollweide projection (a 'flatearth' plot) heatmap. Not used in
manuscript.
- The position along the octiles where the Golgi intersect can be used to generate a heatmap. This also is not in the manuscript
  (except in the methods cartoon).

![windrose and flatearth plots](/Readme_images/Golgi_results_github.png)