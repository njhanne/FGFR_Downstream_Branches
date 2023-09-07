# Downstream branches of FGFR signaling code
Collection of python and R scripts used to analyze data and generate figures in the manuscript

Each directory contains the code needed to analyze a different outcome measure. Each can be used independently from each other, they don't require one another except for the 'DirFileHelpers'.

# Include info here on setting up conda for python...
Much of the code runs in python. I recommend using conda for managing the environment and installing needed packages. It can be downloaded [here]([url](https://www.anaconda.com/download)). When installing if it asks about 'adding python to path' make sure to select 'no'.

Fire up anaconda and create an environment by typing:
```
conda create --name Branches python=3.8
```
When it's done downloading and installing you need to 'activate' the environment by typing:
```
conda activate Branches
```
Now you are in the environment you just created. You can install packages here and they will only be installed on this local environment but not your whole machine. It is much easier to manage your projects this way. You install packages by typing:
```
conda install packagename
```

Before installing any other dependencies I recommend installing Cellpose. The proliferation and Golgi analysis uses Cellpose to segment microscopy images. Detailed explanation of how to install Cellpose for different users is beyond the scope of this littlee readme, but I would recommend visiting their [page here]([url](https://github.com/MouseLand/cellpose)https://github.com/MouseLand/cellpose) and follow their directions to install it locally using conda. Using pip to install packages doesn't always play nice with conda, so I first install their listed dependencies using conda, then install the gui version using pip as they describe. If your computer has an Nvidia gpu then skip down to that section and install cuda pytorch before installing the other dependencies. If conda gets stuck on this, it may be easiest to just start with a clean environment or a dedicated environment for running Cellpose.


# Include info here on setting up Slicer, R, ImageJ
