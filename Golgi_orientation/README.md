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
   of a line drawn between the left and right globular processes on a low magnification 'overview' image of the tissue section.
   See the example image below:
   ![Example of angle adjustment drawn on an overview image](/Golgi_orientation/angle_example.png)