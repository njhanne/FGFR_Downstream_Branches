dir = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/Branches_activation/analysis/raw/";
output_dir = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/Branches_activation/analysis/deconvolved/";
fileList = getFileList(dir);

//activate batch mode
setBatchMode(true);

// LOOP to process the list of files
for (i = 0; i < lengthOf(fileList); i++) {
	// define the "path"
	// by concatenation of dir and the i element of the array fileList
	current_imagePath = dir+fileList[i];
	// check that the currentFile is not a directory
	if (!File.isDirectory(current_imagePath)){
		// open the image and split
		open(current_imagePath);
	currentImage_name = getTitle();

    run("Colour Deconvolution2", "vectors=[H DAB] output=8bit_Transmittance simulated cross hide");
    desiredImage = currentImage_name + "-(Colour_2)";

    selectImage(desiredImage);
    // now we save the DAB image
    currentImage_name = substring(currentImage_name,0,lengthOf(currentImage_name)-4);
    print(currentImage_name);
    saveAs("tiff", output_dir + currentImage_name + "_DAB_deconvolved");
	// make sure to close every images before opening the next one
	run("Close All");
	}
}
setBatchMode(false);
