# Downstream branches of FGFR signaling code
Collection of python and R scripts used to analyze data and generate figures in the manuscript

Each directory contains the code needed to analyze a different outcome measure. Each can be used independently from each other, they don't require one another except for the 'DirFileHelpers'.

# Getting python installed and ready
1. Much of the code runs in python. I recommend using conda for managing the environment and installing needed packages. It can be downloaded [here](https://www.anaconda.com/download). When installing if it asks about 'adding python to path' make sure to select 'no'.

   If you are on Windows, open up anaconda from the start menu. If Macintosh or Linux then just open a terminal. Create an environment by typing:
   ```
   conda create -n Branches python=3.8
   ```
   When it's done downloading and installing you need to 'activate' the environment by typing:
   ```
   conda activate Branches
   ```
   Now you are in the environment you just created. You can install packages here and they will only be installed on this local environment but not your whole machine. It is much easier to manage your projects this way. You install packages by typing:
   ```
   conda install packagename1 packagename2 etc
   ```
2. Before installing any other dependencies I recommend installing Cellpose. Conda tends to get stuck if you don't install this first. Detailed explanation of how to install Cellpose for different users is beyond the scope of this little readme, but I would recommend visiting their [page here](https://github.com/MouseLand/cellpose) and follow the directions to install it locally using conda.
   - Using pip to install packages doesn't always play nice with conda, so I first install their listed dependencies using conda, then install the gui version using pip as they describe.
   - If your computer has an Nvidia gpu then skip down to that section and install cuda pytorch before installing the other dependencies.
   - Here's a brief example that works for me in late 2023 but is likely time and machine dependent! The order of installation is important!
     - Install the dependencies before running installing cellpose with pip. This prevents pip from managing the dependencies.
       ```
       conda install pytorch numpy numba scipy natsort opencv tifffile
       ```
     - Install gui version of cellpose with pip:
       ```
       pip install cellpose[gui]
       ```
     - (Optional, only run if you have an Nvidia gpu). This step can take a while but will save you a lot of time later. Check their directions for more specifics, but visit nvidia's site and download the CUDA toolkit, version 11.6-11.8 I know work. To see which version you have installed, open a terminal window and type:
       ```
       nvcc --version
       ```
       Whatever version shows up there is what you will type into back in the anaconda window. Below I am installing 11.6, but you should match whichever version number is installed on your system:
       ```
       conda install pytorch pytorch-cuda=11.6 -c pytorch -c nvidia
       ```
     - Cellpose should be installed – test it out by typing:
       ```
       cellpose
       ```


# Include info here on setting up Slicer, R, ImageJ
