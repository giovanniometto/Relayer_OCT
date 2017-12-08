# Relayer_OCT

Author: Giovanni Ometto

Originally designed and developed in December 2016.

## Description

Code for the segmentation and analysis of OCT scans

## Dependencies

The following Matlab toolboxes are required:

* Image Processing
* Curve Fitting Toolbox
* Statistics and Machine Learning Toolbox

## Entry points

### Single volume

Read in a folder with image files of the scans:
```
[ILM, RPE, ISOS, THICKNESS] = processVolumeRELAYER(sourceFolder, machineCode, destinationFolder, verbose)
```
Process the volume to extract the ILM, RPE, ISOS and THICKNESS, save csv files in the destination folder with images of the scans and segmentation overimposed.
