Whalen et al. 2020
===========

Code and data to replicate figures found in Whalen et al. 2020, "Delta oscillations are a robust biomarker of dopamine depletion severity and motor dysfunction in awake mice", Journal of Neurophysiology  
https://journals.physiology.org/doi/abs/10.1152/jn.00158.202

Data for figures 1-4, 7 and 9 are included as .mat or .csv (for figure 4 only) files. Other data (with concurrent movement, ECoG or pharmacological interventions) are too large to include in this repository, but can be requested from timcwhalen (AT) gmail.com


Execution
=========

For MATLAB code:  
Ensure the entire Whalen2020 directory and its subdirectories are in your MATLAB path (right click, "Add to Path > Selected Folders and Subfolders"). Then, execute a .m file of the form "Whalen2020_FigureX.m". 

For Python code:  
See Whalen2020_Fig4.py


File and Directory Structure
=========

Files which produce the figures from Whalen et al. are of the form Whalen2020_FigX.m. In general, these load the necessary data from a .mat file in the /data directory into a struct (called "data") and call a file with the suffix "\_batch", then plot the results. These batch scripts take the "data" struct and add a substruct containing the results and parameters used (for ease of replication). You can specify your own parameters by adding the substruct and setting the parameters yourself, as explained in each batch script's notes.

All other functions in the main directoy serve more general purposes and may be useful beyond the scope of this paper (but see "Looking for more?" below for a more useful repository containing this code)

/helpers contains specific helpers for plotting and managing data in batch scripts and are probably not very useful more generally.

/ext_depend contains external dependencies which were obtained from MATLAB Exchange or directly from colleagues. They were not written by the authros of this paper; their attributions are in the code.


Looking for more?
=========

This repo is primarily for replicating the figures in Whalen et al. 2020, and thus it is not updated to freeze these files in time. For a more up-to-date repo which includes this code (except those files prepended Whalen2020) and other useful code for analysis of neural data, see https://github.com/tcwhalen/InVivo

Special thanks to Dr. Rob Turner for several dependencies and other code on which some of these functions are based.