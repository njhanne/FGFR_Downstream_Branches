# 3D geometric morphometrics
This code is used to perform shape analysis on 3D reconstructions of uCT scanned embryos.
The code has some scripts that help with converting the files from Scanco aim files to ply volumes, too.

## Preparing the samples for landmark placement
As mentioned on the main page - the landmark placement is done in [3DSlicer](https://download.slicer.org/), a free 3D volume program.
I recommend also getting the [SlicerMorph Project](https://slicermorph.github.io/) addons which help with converting the Scanco volumes.

### Opening Scanco .aim quickly in Slicer
Inside the 'Slicer_Scripts' folder the 'AIM_to_NHDR.py' script will create a header file that can be used to open .aim files in Slicer.
There was a blogpost by Murat Maga detailing how to generate these nrrd header files, but it's been deleted. I've copied it here for convenience:

```
Importing microCT data from Scanco into Slicer

In this tutorial, I will show how to import a ‘raw image’ data into Slicer directly using Nearly Raw Raster Data (NRRD)
format specification. I will be using a microCT data sample from a Scanco scanner courtesy of Dr. Chris Percival,
however steps are equally applicable to raw image data from any other scanner vendor. 
For this to work, we need to know the dimension, data type and other metadata associated with the imaging session. 
These are typically provided as a log file. All files associated with this tutorial can be found in here.

Log file contains all sorts of potentially useful information about the scanning session, but the things we are 
interested are located in the very first section: Where image data starts in the raw volume (byte offset),
XYZ dimensions (622 642 1085), data type (short, or 16 bit), and voxel size (0.02, 0.02, 0.02). 
These are key pieces of information that is necessary to describe the image data correctly through a ‘detach header’ 
(for more information see http://teem.sourceforge.net/nrrd/format.html).

Once you download the sample .aim file, open a text editor and paste the following lines, and save it in the same folder 
as the file you download with extension .nhdr (e.g, test.nhdr).

NRRD0004
type: short
dimension: 3
sizes: 622 642 1085
spacings: 0.02 0.02 0.02
byteskip: 3195
data file: ./PERCIVAL_370.AIM
endian: big
encoding: raw

If you drag and drop the nhdr file into Slicer, you should see the sample data being loaded (this is a fairly large 
dataset, may take sometime). If this is the first time you will be working with this dataset, you are done. However, 
if you have additional data (e.g., landmarks) collected outside of Slicer, it is important to check whether the coordinate 
systems are consistent between two different programs. Dr. Percival also provided a set of landmarks coordinates that was 
collected from this specimen using Amira. If you drag and drop the landmarks.fcsv into Slicer, you will see that landmarks 
do not correctly line up.

We need to identify the coordinate system the original landmarks were collected, and describe in the NHDR file. This 
step will involve some troubleshooting (e.g., collecting a few of the landmarks in Slicer, and comparing reported 
coordinates.) In this case, all positive values of the original landmark coordinates suggested to me that it is likely 
the original coordinate system was Right-Anterior-Superior (or RAS in short). I modified the header file to include this 
as the space tag. One other difference is that you need to remove the spacings tag, and then redefine it as the space directions.

NRRD0004
type: short
dimension: 3
sizes: 622 642 1085
byteskip: 3195
data file: ./PERCIVAL_370.AIM
endian: big
encoding: raw
space: RAS
space directions: (0.02,0,0) (0,0.02,0) (0,0,0.02)

If you change the content of the NHDR file with the code above, and reload your volume and render, now you should see all 
the landmarks line up correctly. It is always good to familiarize yourself with the coordinate systems of the digitization 
program you are using. You can find a discussion of coordinate systems here.


This entry was posted in 3D Imaging, data import, Slicer on September 6, 2018 by maga.
```

In our case we have some info needed in the .aim file, and some other info in the .hx file. I think the hx is for importing into
Amira? The .text from ipl (image processing log, maybe?) file is generated by Scanco to help process the larger scan into individual samples - it also has some useful
information needed in the header. Unfortunately, there is still more information required not found in these files to generate a working .nrrd header file.
We can get the remaining info using the [Raw Image Guess](https://github.com/acetylsalicyl/SlicerRawImageGuess) plugin.

In any case, the python script I wrote should generate good .nhdr files for the scans used in this project. You can open the .nhdr file in Slicer
and it will open the corresponding .aim file!

### Segment .aim and write .ply
Inside the 'Slicer_Scripts' folder the 'load_and_segment.py' script will load each of the .aim files in a directory into Slicer 
(sequentially, not all at once), run some segmentation steps, and write the segmentation as a .ply file. The .ply files are much smaller than the .aim files,
making them faster to work with in following steps. 

The script is run inside of Slicer but is launched from the system terminal. It uses its own python environment, not your conda env. 
Luckily the code here doesn't need any special libraries so it will work fine.

The bad news is that you need to change directory in your terminal to the location of the 'load_and_segment.py' file.

On Windows:
``` 
& '/full/path/toSlicer.exe' --python-script "./load_and_segment.py" --no-splash
```
On Macintosh / Linux:
```
/Applications/Slicer.app/Contents/MacOS/Slicer --no-splash --python-script "/full/path/to/FGFR_Downstream_Branches/Geometric_morphometrics/3D/Slicer_Scripts/load_and_segment.py"
```

## Placing landmarks
I first place the 12 fiducial landmarks in the order shown in the figure below. Next I make a semi-curve and copy-paste in the 3rd and 4th fiducial landmark. I then move the 9th fiducial landmark so that it lies on the line drawn between the 3rd and 4th.

![fiducial_landmarks](/Readme_images/Fiducial_layout.png)

Next I make another semi-curve and copy-paste in the 9th and 10th fiducial landmark. In the curve settings dropdown I select 'shortest distance on surface' and constrain the model to the current volume. Then I resample the line to 5 evenly spaced landmarks.

![center_resampling](/Readme_images/center_semicurve.jpg)

Next I draw 4-5 points around the outside of the nasal pit and resample the line to 9 evenly spaced landmarks. Repeat for the other side of the face.

The semi-landmarks are a grid of 3x3 landmarks from the nasal pit to the midline of the FNP. The 'rows' follow a similar spacing as the midline curve-semilandmarks. Repeat for the other side of the face. The below figure shows all fiducial landmarks, the center curve semilandmarks, and the right side semilandmarks.

![right-side-lms](/Readme_images/semi-landmarks.png)

## Analysis
Almost the entirety of the analysis code was adapted from Marta Vidal Garcia. Git blame will show 'Nicholas Hanne'
as the editor for many lines, but they are mostly edits to the existing work of Marta. Thank you! All of the analysis
is performed in R using [Morpho](https://cran.r-project.org/web/packages/Morpho/Morpho.pdf) and [geomorph](https://cran.r-project.org/web/packages/geomorph/geomorph.pdf).
Explaining GM analysis is maybe beyond the scope of this readme, but our R scripts are well commented (imo) and the 
official documentation linked in the previous sentence are extremely helpful.

'GMM_analyses.R' should be run before 'Asymmetry.R'.

![GMM_analysis.R output](/Readme_images/3D_GM.png)