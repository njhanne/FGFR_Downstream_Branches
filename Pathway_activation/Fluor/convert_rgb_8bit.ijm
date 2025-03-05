// https://gist.github.com/romainGuiet/cf42f3b1d31222a76d602dfe2f028894
dir = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/Branches_activation/original_images/Fluor/";
output_dir = "C:/Users/njhan/Box/FGF_inhibitor_paper_5-26-2020/data/Branches_activation/analysis/raw/Fluor/";
fileList = getFileList(dir);

//activate batch mode
setBatchMode(true);

// LOOP to process the list of files
for (i = 0; i < lengthOf(fileList); i++) {
	// define the "path" 
	// by concatenation of dir and the i element of the array fileList
	current_imagePath = dir+fileList[i];
	// check that the currentFile is not a directory
	if ( !File.isDirectory(current_imagePath)) {
		run("Bio-Formats Macro Extensions");
		Ext.openImagePlus(current_imagePath);
		
		currentImage_name = getTitle();
		print(currentImage_name);
		currentImage_name = replace(currentImage_name, "(.+\/)(.+)\.", "$2");
		currentImage_name = substring(currentImage_name,0,lengthOf(currentImage_name)-3);
		print(currentImage_name);

		if (bitDepth() != 24) {
			run("RGB Color");	
		}
		
		run("8-bit");
		saveAs("tiff", output_dir + currentImage_name); 
		run("Close All");
	}
}
setBatchMode(false);