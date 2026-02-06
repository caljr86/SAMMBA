// This macro is made of loops to screen several .JPG images from commercial grade cameras.
// This assumes that cameras cannot name the files diretly according to the well positions

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
experiment = "breeding/";

//Input folder
raw = root + experiment + plate + "/" + date + "/raw/";
// Output cropped raw images
cropped = root + experiment + plate + "/" + date + "/cropped/";
//Output path for segmented images
outputDir = root + experiment + plate + "/" + date + "/segmented/";
//Input path for LabKit model
labkitModelDir = root + "script/"+ model;
//Input path for plate_design.csv file
design_path = root + experiment + plate + "/plate_design.csv";

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
  				fileName = list1[i];
  
   	 			// Open each image and duplicate for croping and scaling
    				open(raw + fileName);
    				img = getTitle();
    				run("Duplicate...", " ");
    					//Set scale before croping
							run("Set Scale...", "distance=5206.5160 known=2700 unit=um");
    					//setTool("rectangle");
							makeRectangle(978, 18, 4440, 4142);
								run("Crop");
							saveAs("ZIP", cropped + fileName);
  								close("*");
  
				// Segmenting and analyzing with LabKit 

	//Create objects to change file names when opening
   img0 = IJ.pad(i+1, 4);
   img1 = "IMG_" + img0 + ".zip";
  
	//Open cropped images and remove the extensions
   open(cropped + img1);
    img = replace(img1, ".tif", "");
    img = replace(img1, ".zip", "");
    
    //Get the well name from list according to image sequence
    line = wells[i+1];
	line2 = split(line, ",");
	well = line2[2];

	//Create path to modify the name of the segmented image and results, including the right well according to design file

segmentedPath_zip = outputDir + "zip/"+ date + "_" + img + "_" + well + "_segmented";
segmentedPath_res = outputDir + "results/" + date + "_" + img + "_" + well + "_segmented";
    
				// apply labkit model
   					run("Segment Image With Labkit", "segmenter_file=" + labkitModelDir + " use_gpu=false");
   					saveAs("ZIP", segmentedPath_zip);
					close();
close();	

//Open segmented image and analyze
  open(segmentedPath_zip + ".zip");

wait(1000);

				// Threshold image for Analyze Particle function
					setAutoThreshold("Default dark");
					//run("Threshold...");
					//setThreshold(0, 0);
						setOption("BlackBackground", true);
						run("Convert to Mask");
						wait(1000);

   				// Get count and area of particles using 5um as lower threshold
//Change this threshold based on previous analysis.
					run("Analyze Particles...", "size=5-Infinity display");
   
					//save as zip to compress and save metadata

						saveAs("Results", segmentedPath_res + "_Results.csv");
   
    print("Saved Image and Results: " + segmentedPath);
    print("i = " + i);

//close all windows
close("*");
close("Results");
 //Dump memory usage for next image
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





















