// This macro is made of loops to screen several .vsi images from Olympus Microscope
// This assumes that software can name the files diretly according to well positions

// Record the start time
startTime = getTime();

//Type images date
model = "plate1.classifier";

//Add inside the brackets the folder(s) you want to analyse
plates = newArray("plate1","plate2","plate3")

//Loop to screen through desired plates
for(p= 0; p < lengthOf(plates); p++) {
	plate = plates[p];
print("Processing: " + plate);

// Loop to analyze several folders in a row per plate. Add here the date(s) you want to analyse
folders = newArray("260109", "260112", "260114", "260116", "260119", "260121", "260123");

//Loop to screen through desired dates inside each plate
for(f = 0; f < lengthOf(folders); f++) {
	date = folders[f];
print("Processing folder: " + date);

// ImageJ macro for croping, scaling and segmenting a list of images
// Check or reorder the folders according to your preferred organization

root = "/media/Disk1/user1/";
//Input folder
raw = root + plate + "/" + date + "/raw/";
// Output cropped raw images
cropped = root + plate + "/" + date + "/cropped/";
//Output path for segmented images
outputDir = root + plate + "/" + date + "/segmented/";
//Input path for LabKit model
labkitModelDir = root + "script/"+ model;
//Input path for plate_design.csv file
design_path = root + plate + "/plate_design.csv";

//Open plate design file
design = File.openAsString(design_path);

//Extract from design object an array separated by lines, which contain the well names
wells = split(design, "\n");

//Get list of all files and folders in raw path
list = getFileList(raw);


// Record the start time
startTime = getTime();
 
 //Loop to segmente images individually

		// In this part of the loop you can change the "0" and "list.length" by any number according to the total number of files
			for (i = 0; i < list.length; i++) {
  				fileName = list[i];

run("Bio-Formats Importer", "open="+ raw + fileName + " color_mode=Default rois_import=[ROI manager]"+
" view=Hyperstack stack_order=XYCZT series_1");
//run("Z Project...", "projection=[Sum Slices]");
}

img = getTitle();
img1 = replace(img, ".vsi", "");
//img2 = replace(img1, "SUM_img_", "IMG_");
img2 = replace(img1, " - PI", "");

makeRectangle(697, 0, 2691, 2278);
run("Crop");
saveAs("ZIP", croped + fileName);


segmentedPath_zip = outputDir + "zip/"+ date + "_" + img2 + "_segmented";
segmentedPath_res = outputDir + "results/" + date + "_" + img2 + "_segmented";

 // apply labkit model
   run("Segment Image With Labkit", "segmenter_file=" + labkitModelDir + " use_gpu=false");
   saveAs("ZIP", segmentedPath_zip);
close();
close();

//Open segmented image and analyze
  open(segmentedPath_zip + ".zip");

wait(1000);

setAutoThreshold("Default dark");
//run("Threshold...");
//setThreshold(0, 252);
setOption("BlackBackground", true);
run("Convert to Mask");

wait(1000);

   // Get count and area of particles using 5um as lower threshold
run("Analyze Particles...", "size=10-Infinity display");
   
   //Set name and path for saving

//save as zip to compress and save metadata

saveAs("Results", segmentedPath_res + "_Results.csv");
   
    print("Saved Segmented Image and Results Table");
    print("i = " + i);

//close all windows
close("*");
close("Results");
 
   run("Collect Garbage");

// Warning that loop is over
print("################____NEXT LOOP_____################"); 

}
} 
}

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

//Save log for control and debuging
logpath = outputDir + "Log" + i + ".txt";
selectWindow("Log");
saveAs("Text", logpath);





















