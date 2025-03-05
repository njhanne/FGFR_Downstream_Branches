dir = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/Branches_activation/analysis/cellpose/Fluor/expanded_mask/";
output_dir = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/Branches_activation/analysis/raw/Fluor/activation/";
fileList = getFileList(dir);

//activate batch mode
setBatchMode(true);

// LOOP to process the list of files
for (i = 0; i < lengthOf(fileList); i++) {
	current_imagePath = dir+fileList[i];
	print(current_imagePath);
	// check that the currentFile is not a directory
	if (!File.isDirectory(current_imagePath)){
		// open the image and make the rois
		open(current_imagePath);
		currentImage_name = getTitle();
		run("Label Map to ROIs", "connectivity=C4 vertex_location=Corners name_pattern=r%03d");

		roiManager("Select", 1);
		run("Select All");

		currentImage_name = substring(currentImage_name,0,lengthOf(currentImage_name)-4);
		print(currentImage_name);
		roiManager("Save", output_dir + currentImage_name + ".zip");

		// make sure to close every images befores opening the next one
		roiManager("Select", 1);
		run("Select All");
		roiManager("Deselect");
		roiManager("Delete");
		run("Close All");
	}
}
setBatchMode(false);