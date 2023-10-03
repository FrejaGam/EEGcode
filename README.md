# EEGcode
Repository for matlab code from [Kaslab](https://research.rug.nl/en/organisations/kas-lab-behavioural-neuroscience) at Rijksuniversiteit Groningen. 
The code here is designed to be used on .EDF+ files in general, but adapted to files recorded in TainiLive.

# User Guide
Below is the description of the script for single ERPs along with some general features and dependencies for the data analysis of local field potentials (LFP). For recordings of spontaneous activity the preprocessing and filters are the same as for ERPs thereafter the pipeline forks [scoring](https://github.com/FrejaGam/EEGcode/wiki/Scoring-of-vigilance-states). More data processing scripts can be found in the [wiki](https://github.com/FrejaGam/EEGcode/wiki).

## EEG processing workflow
This document describes the general workflows, with suggested software and scripts for processing EEG recorded either with stimulus (event related potentials) or without stimulus (resting state, sleep).
For all files the processing is divided into three parts
1.	preprocessing is done separately for each animal.
2.	Analysis and average within groups is done in a script separate from preprocessing.
3.	Statistical comparisons between groups are typically done in R.

Event related potentials (incl. VEP, AEP, gating, SSVEP and chirp)
For all event-driven paradigms the scripts contain averaging over events.
The matlab scripts for these paradigms require the toolboxes: communications, DSP system, image processing and signal processing.
They further require the matlab functions: 
- blockEDFload.m [https://github.com/DennisDean/]
- parse_taini_edf_notations.m [https://tainitec.com/]
- filterEEG.m *(this function applies both low-pass and high-pass Butterworth filters)*

### VEP and AEP (ERP)
*Matlab (preprocessing)*
- Open the script 1_1_fullfile_SingleERP.m
- The first block of the script asks where the data file is and loads it to matlab
- The second block makes a figure of the triggers in the data and numbers them *(the trigger numbers don’t vary between trials for VEP, but they might for AEPs, it’s something with Avisoft )*
This is used if there are multiple stimulus-paradigms in the data file. It is then possible to only choose the triggers relevant for this particular analysis.
- The third block is where the EEG is cut around event points chosen in the second block. The specific recording channels are also added here. A matrix is made with the channels in a specific time interval around the trigger.
- In the fourth block the data is filtered, artifact rejected and averaged. A matrix is created, containing the number of artifact rejected trials for each channel.
- The fifth block is a figure of the raw data and the filtered data. If the filtered data is a straight line, then the thresholds for artifact rejection are most likely set wrong.
- The sixth block is a figure of the ERP for each channel.
- The last block is for saving the ERP data and printing the number of artifact rejected trials (can be copy-pasted to an excel file). The ID variable is set to extract the animal number from the name of the .edf file. The path1 variable can be used to specify the path for the saved file. The variable IDs can be used to add a suffix related to experiment such as batch or other relevant info. The suffix should be the same for data that needs to be averaged.

*Matlab (group level analysis)*
- Open the script 1_2_Averaging_ERP.m (this script calls the function [shadederrorbars.m](https://nl.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar))
- The first block is for loading the data relevant for the averaging these should have been given the same suffix. The sampling frequency and the relevant time interval are hardcoded in this block.
- In the second the block the data channels of the different files, are divided into relevant genotypes or treatment groups. And figures are made of the average ERP within each group.
- For ERP data the amplitude and time should be extracted separately. In terms of matlab programming this is done be extracting the max. and min. values of relevant time windows. These are defined in the third block of the script. This loop will draw the cursors, delimiting the time windows, onto each channel in a specific group. This block only draws the cursors the extraction is done in the following block.
- The fourth block the amplitude and index are extracted for each peak. The amplitudes go into the matrix peaksA. The indices go to peaksI. 
These indices have to be added to the cursor to get the data point e.g., if the index is 15 and the cursor is at 80 then the actual index for that peak is happening at the 95th data point in the time interval. To get the time this value should be divided by the sampling rate.
The data can be copy-pasted from peaksA and peaksI to a spreadsheet (e.g., excel). From excel it can then be loaded into a statistical program of choice. 
Remember experiment relevant variables such as ID, batch, side, electrode etc.

