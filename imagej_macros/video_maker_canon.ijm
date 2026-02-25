// Record the start time
startTime = getTime();

plates = newArray("plate58");
well = "C21";

for(p = 0; p < lengthOf(plates); p++) {
	plate = plates[p];
print("Processing: " + plate);

root = "/home/cicero/Dropbox/1_posdoc/1_CCMAR/experiments/384WP/biobank/growthrate/";
timepoints = root + plate + "/timepoints/";

// Loop to analyze several folders in a row
//folders = newArray("260127","260129","260202","260204", "260206");
folders = getFileList(timepoints);

for(f = 0; f < lengthOf(folders); f++) {
	date = folders[f];
print("Processing folder: " + date);

// ImageJ macro for croping, scaling and segmenting a list of images

raw = timepoints + date + "/raw/";

outputDir = root + plate;

//Input path for plate_design.csv file
design_path = outputDir + "/plate_design.csv";
//Open plate design file
design = File.openAsString(design_path);
//Extract from design object an array separated by lines, which contain the well names
wells = split(design, "\n");

//Loop to get only .vsi. Choose the right i based on your image sequence

for (i = 0; i < wells.length; i++) {
	 if (indexOf(wells[i], well) != -1) {
        foundIndex = i;
        break; // stop at first match
    }
}

wellline = wells[foundIndex];
wellline2 = split(wellline, ",");
wellnumber = wellline2[3];
wellnumber2 = parseInt(wellnumber);

img0 = IJ.pad(wellnumber2 + 1, 4);
filename = "IMG_" + img0 + ".JPG";

open(raw + filename);
run("Brightness/Contrast...");
	setMinAndMax(0, 60);
}}

run("Images to Stack", "use");

	run("Linear Stack Alignment with SIFT", "initial_gaussian_blur=1.60 steps_per_scale_octave=3 minimum_image_size=64 " +
	"maximum_image_size=1024 feature_descriptor_size=4 feature_descriptor_orientation_bins=8 closest/next_closest_ratio=0.72 " +
	"maximal_alignment_error=1 inlier_ratio=0.1 expected_transformation=Rigid interpolate");

StackName = plate + "_" + well + "_" + img0;
selectWindow("Stack");
run("AVI... ", "compression=JPEG frame=1 save=" + outputDir + "/"+ StackName);
close("*");


// Record the end time
endTime = getTime();

// Warning that processing is over
print("################____ACABOU!!!_____################"); 

 // Happy gametophyte drawing
print("_____________________________________░___________________________________________________");
print("_____________________________________░░___________░░_______________________________");
print("________░________░___________░░____░░░________░░_________________________________");
print("__________░░______░░░_______░░____░░░_____░░░_____░░__________________________");
print("____________░░_______░░__░░░░__░░░░__░░░░_____░░___________________________");
print("_____________░░░_______░░░░░░░░░░░░░░░░░░░░____░░____________________________");
print("_________________░░░░__░░░░░░░░░░░░░░░░░░_____░░_____________________________");
print("________________░░░░░░░░░░░░░░░░░░░░░░░░░░░░______________________________");
print("_____________________░░░░░░░░░░░░░░░░░░░░░░__________________________________");
print("____________░░░░░░░░░░░░██░░░░██░░░░░░░░░░░░░░__________________________");
print("_____________________░░░░░░░░░░░░░░░░░░░░░░░░░░______________________________");
print("__________________░░░░░░░░░░░░░░░░░░░░░░░░░░________________________________");
print("______░░░░░░░░░░░░░░░░██░░░░██░░░░░░░░░░░░░░__________________________");
print("___________________░░░░░░░░░█████░░░░░░░░__________________________________");
print("__________________░░░░░░░░░░░░░░░░░░░░░░░░________________________________");
print("___________________░░░░░░░░░░░░░░░░░░░░░░___________________________________");
print("________________░░______░░░░░░░░░░░░░░_____░░_____________________________");
print("_______________░░_______░░____░░░░__░░___________░░____________________________");
print("_____________░░________░░______░░_______░░__________░░________________");
print("____________________░░░________░░__________░░________░░____________________");
print("_________________░░░░__________░░____________░░░_________________________") ;


// Calculate and print the elapsed time
elapsedTime = (endTime - startTime) / 1000;
print("Batch processing completed in " + elapsedTime + " seconds.");
print("Batch processing completed in " + elapsedTime/60 + " minutes.");
print("Batch processing completed in " + elapsedTime/3600 + " hours.");





















