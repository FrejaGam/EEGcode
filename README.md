# EEGcode
Test repository for matlab code (some for R) from Kaslab. 
The code here is designed to be used on .EDF+ files created with a TainiTec system.

# User Guide
## EEG processing workflow
This document describes the general workflows, with suggested software and scripts for processing EEG recorded either with stimulus (event related potentials) or without stimulus (resting state, sleep).
For all files the processing is divided into three parts
1.	preprocessing is done separately for each animal.
2.	Analysis and average within groups is done in a script separate from preprocessing.
3.	Statistical comparisons between groups are typically done in R.

Event related potentials (incl. VEP, AEP, gating, SSVEP and chirp)
For all event-driven paradigms the scripts contain averaging over events.
The matlab scripts for these paradigms require the toolboxes: communications, DSP system, image processing and signal processing.
They further require the matlab script: blockEDFload.m, parse_taini_edf_notations.m, filterEEG.m

### VEP and AEP (ERP)
*Matlab (preprocessing)*
- Open the script fullfile_SingleERP.m
- The first block of the script asks where the data file is and loads it to matlab
- The second block makes a figure of the triggers in the data and numbers them (the trigger numbers don’t vary between trials for VEP, but they might for AEPs, it’s something with Avisoft )
This is used if there are multiple stimulus-paradigms in the data file. It is then possible to only choose the triggers relevant for this particular analysis.
- The third block is where the EEG is cut around event points chosen in the second block. The specific recording channels are also added here. A matrix is made with the channels in a specific time interval around the trigger.
- In the fourth block the data is filtered, artifact rejected and averaged. A matrix is created, containing the number of artifact rejected trials for each channel.
- The fifth block is a figure of the raw data and the filtered data. If the filtered data is a straight line, then the thresholds for artifact rejection are most likely set wrong.
- The sixth block is a figure of the ERP for each channel.
- The last block is for saving the ERP data and printing the number of artifact rejected trials (can be copy-pasted to an excel file). The ID variable is set to extract the animal number from the name of the .edf file. The path1 variable can be used to specify the path for the saved file. The variable IDs can be used to add a suffix related to experiment such as batch or other relevant info. The suffix should be the same for data that needs to be averaged.
