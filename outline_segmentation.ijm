// Florian Stroehl    fs417@cam.ac.uk
// 2016-10-21
// v1.0

// get name
name = getTitle;
dir = getDirectory("image");

dotIndex = indexOf(name, ".");
title = substring(name, 0, dotIndex); 

// manual threshold
run("Threshold...");
waitForUser("set the threshold and press OK, or cancel to exit macro"); 
run("Convert to Mask");

path = dir + title + "_MASK.tif"; 
saveAs("Tif", path);
close()
