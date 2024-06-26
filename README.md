# mycroplanter

This repository contains the processed data and code for the manuscript, "Distinguishing deterministic and stochastic processes in the rhizosphere microbiome" by Chen et al (anticipated publication 2024). 

Raw data (image files, plate maps, etc) can be found on Dryad (DOI: 10.5061/dryad.w9ghx3fxd). Only processed image data (pixel data) and cleaned metadata files are found here. The data are split into four analysis paths:


* Competition ratio experiment (compratio)
* Fluorescence experiment (fluor)
* Priority effects experiment (priorityeffets)
* Priority effects on plates (priorityplates)

In the main directory, you will find the image processing script, 'process_athal_scan.py'. IMPORTANT NOTE: For the manuscript, we use 'process_athal_scan.py' in the main directory of this repositiory, which is a simplified version of the full script. THe full (and updated) version of the script is in "mycroplanter_processing_tutorial". Both scripts are both capable of the same fucntions, but the latter has more options for the user and a streamlined "default" pipeline. The latter also uses a custom progress bar instead of the pre-built package 'ProgressBar'. If you are using this program for your own research, please use the version in the tutorial.





