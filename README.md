# SMTI
Software to analyse data generated using single molecule translation imaging.

3rd party programs required:
 - rapidSTORM (http://www.super-resolution.biozentrum.uni-wuerzburg.de/research_topics/rapidstorm/)
 - Fiji/ImageJ (https://fiji.sc/)
 - Matlab
Step 1: Localise SMTI data sets with rapidSTORM using the rapidSTORM configuration file "rapidSTORM_configuration_file.txt".
        Dependent on the used microscope the parameters in the file might have to be adapted.
Step 2: Create an outline image of a fluorescence widefield image of the sample using the macro "outline_segmentation.ijm" in Fiji.
Step 3: Open "SMTI_Analysis.m" in Matlab and set the input system parametes according to your system. 
        Setting the respective flags will create time courses, rate histograms, density maps, localisation maps, measure rates, 
        and create movies from your localisation data file (produced by rapidSTORM) and your outline image (produced by Fiji)
