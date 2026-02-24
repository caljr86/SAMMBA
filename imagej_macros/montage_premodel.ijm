// Record the start time
startTime = getTime();

plates = newArray("plate78")

for(p = 0; p < lengthOf(plates); p++) {
	plate = plates[p];
print("Processing: " + plate);

// Loop to analyze several folders in a row
//folders = newArray("260127","260129","260202","260204", "260206");
folders = newArray("260211", "260213", "260217", "260218", "260220");

for(f = 0; f < lengthOf(folders); f++) {
	date = folders[f];
print("Processing folder: " + date);

// ImageJ macro for croping, scaling and segmenting a list of images

root = "/home/cicero/Dropbox/1_posdoc/1_CCMAR/experiments/384WP/aquadiversify/growthrate/";
raw = root + plate + "/timepoints/" + date + "/raw/";

outputDir = root + plate;

list = getFileList(raw);

//Loop to get only .vsi. Get only H11, H12 and H13

for (i = 0; i < list.length; i++) {
	if(!endsWith(list[i], ".vsi"))
	continue{

fileName = list[i];

run("Bio-Formats Importer", "open="+ raw + fileName + " color_mode=Default rois_import=[ROI manager]"+
" view=Hyperstack stack_order=XYCZT series_1");
//run("Z Project...", "projection=[Sum Slices]");
}}}

run("Images to Stack", "use");
run("Make Montage...", "columns=" + list.length/2 + " rows=" + folders.length + " scale=1");
//run("Brightness/Contrast...");
	//setMinAndMax(0, 20);

StackName = plate + "_" + "Montage";
selectWindow("Montage");
saveAs("ZIP", outputDir + "/" + StackName);
close("*");
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





















