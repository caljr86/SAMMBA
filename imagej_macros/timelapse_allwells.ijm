plate = "plate56";

//Input folder

root = "/home/cicero/Dropbox/1_posdoc/1_CCMAR/experiments/384WP/growthrate/";
outputDir = root + plate + "/stacked_aligned/";
	
//Input path for plate_design.csv file
design_path = root + plate + "/plate_design.csv";
//Open plate design file
design = File.openAsString(design_path);

//Extract from design object an array separated by lines, which contain the well names
wells = split(design, "\n");

 // Loop to analyze several folders in a row
folders = newArray("250623","250625", "250627", "250630", "250702", "250704", "250710", "250721");

for(w = 0; w < wells.length; w++) {
		
				img0 = IJ.pad(w + 1, 4);
				filename = "IMG_" + img0;
				filename1 = filename + ".zip";
				days = lengthOf(folders);
		for(i = 0; i < days ; i++) {
	date = folders[i];
	croped = root + plate + "/" + date + "/croped/";
		
				open(croped + filename1);
				run("Scale...", "x=0.1 y=0.1 width=444 height=414 interpolation=Bilinear average create");
				close(filename + ".tif");
}

run("Images to Stack", "use");
selectImage("Stack");

run("Linear Stack Alignment with SIFT", "initial_gaussian_blur=1.60 steps_per_scale_octave=3 minimum_image_size=64 " +
	"maximum_image_size=1024 feature_descriptor_size=4 feature_descriptor_orientation_bins=8 closest/next_closest_ratio=0.92 " +
	"maximal_alignment_error=25 inlier_ratio=0.05 expected_transformation=Rigid interpolate");

    line = wells[w + 1];
	line2 = split(line, ",");
	well = line2[2];


aligned = "Aligned " + days + " of " + days;
selectImage(aligned);
saveAs("ZIP", outputDir + well + "_aligned.jpg");
close("*");
}