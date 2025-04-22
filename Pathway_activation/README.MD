# Pathway activation
This code quantifies the amount of staining / labeling of DAB or fluorescent IHC images, respectively. The DAB folder has the 
code for quantifying by DAB, which fundamentally does not work and should not be used. But there are some interesting bits of 
code, like the plugin for ImageJ that performs color deconvolution. The Fluor folder is what is used in the manuscript and
is more robust method for this quantification.

# Fluorescent quantification of pathway activation
## Preparing images and segmenting them with Cellpose
These labeled proteins reside primarily outside the nucleus so the technique used to quantify apoptosis and proliferation
will not work here. Instead our strategy will be to mask the nuclei, expand the masks until they abut one another to approximate
the cells' membrane, and then get the average label intensity in the expanded area.

1. The images were captured using the Leica fluorescent microscope in our lab. It is not confocal, the images were all captured
with the 40x objective. The fluorescent images are RGB since they are captured with a camera, not by the scope itself.
Luckily this set of images was captured with the channels completely separated, so converting them to 8bit and forming a 
stack is simple. The imageJ script 'convert_rgb_8bit.ijm' will perform that.

2. I used Cellpose 3 to segment the nuclei. In this case the new 'cyto3' model can segment the nuclei very well and hand
training is not necessary. The python script 'run_cellpose_pathway_activation.py' will generate segmentations of the nuclei.

3. The nuclei labels can be expanded to approximate the edge of the cell membrane. See
[this great resource](https://bioimagebook.github.io/chapters/2-processing/6-transforms/transforms.html#partitioning-images-with-voronoi) here that explains
that explains how to use a watershed algorithm to expand the nuclei masks so they don't overlap. I expanded the masks by 
15 pixels, which is about 1/2 the diameter of the nuclei segmentations. 'nuc_mask_expansion.py' performs this expansion
and saves the new label images.

The below images show 8bit image of the DAPI stained image on the left, the cellpose segmentation of the nuclei in the middle,
and the expanded nuclei segmentations on the right.
![cellpose output](/Readme_images/activation_cellpose_expansion.png)

## Thresholding activated cells
The threshold for positive vs negative cells is different in each image and slightly subjective. The strategy for reducing 
subjectivity was to convert the segmentations into ImageJ rois. These rois can then be selected in ImageJ to measure
the average brightness of the phosphorylated target antibody. I selected 3 cells that I thought were barely positive and
3 cells that were nearly positive, then took their average to determine the threshold brightness for each image.

4. The 'labels_to_roi.ijm' ImageJ script will convert the expanded label images into ImageJ rois.
5. These rois can can then be imported into ImageJ with the stack created from step 1. The second channel in the stack
is the fluorescently labeled phosphorylated downstream targets. Selecting a cell's roi and then clicking 'measure' in ImageJ
will calculate the mean grey value of that cell. I measured 6 cells' brightness in each image and put them in a csv.

The below left image shows the 8bit fluorescence image of one of the phosphorylated target antibodies. The middle image shows
the overlaid ImageJ rois (as created in step 4 above). The right image is a zoomed-in portion of the middle image, with 6
cells highlighted. The red circled cells are 'barely negative' and the green are 'barely positive'. The average brightness
of these 6 cells is used to calculate the threshold for each image.
![threscholding_imageJ](/Readme_images/activation_threshold_pic.png)

## Filtering by distance from frontal ectoderm and quantification
Another consideration is that most of the positive signaling is adjacent to the frontal ectoderm. If we sample the full image
then the vast majority of negative cells farther from the ectoderm will dilute the ratio of positive vs negative cells, 
which could likely mask the impact of the inhibitor treatment. To only look at our region of interest I draw a freehand
line roi in ImageJ and use python to only select cells 200 um away from this line. This also filters out the cells in
the ectoderm from being counted.

6. In ImageJ I draw a freehand line along the boundary between the frontal ectoderm and the mesenchyme. See the left image, below.

7. This roi, the expanded label image, and the 8bit image of fluorescent brightness are analyzed in 'mask_counting.py'.
As demonstrated in the above rightmost image, the threshold is made by averaging the brightness in the selected barely
positive and barely negative cells. This threshold is then applied to calculate a ratio of positive to negative cells in 
the entire image. Finally, the ImageJ freehand line roi is converted into a shape object, which is then extended 200 um
into the mesenchyme, further filtering out cells too deep in the tissue section. These ratios are saved as a .csv file.The code also saves a diagram showing the filtered cells. An example is shown on the right image, below.

8. Finally, the csv is loaded and analyzed in R with the 'phosphorylation_quant.R' script. This script is straightforward
and annotated.

![threscholding_roi_distance](/Readme_images/activation_threshold_distance.png)