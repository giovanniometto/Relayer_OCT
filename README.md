# Relayer_OCT

Author: Giovanni Ometto

Work address: C274 Tait Building City, University of London, London, EC1V 0HB (UK)

Email: giovanni.ometto@city.ac.uk

Website: http://www.city.ac.uk

## Description

Code for the segmentation and analysis of OCT scans

## Dependencies

The following Matlab toolboxes are required:

* Image Processing
* Curve Fitting Toolbox
* Statistics and Machine Learning Toolbox

## Entry points

### Single volume

Read in an OCT volume from HDF5 file:
```
[ILM, RPE, ISOS, THICKNESS] = processVolumeRELAYER(sourceFolder, machineCode, destinationFolder, verbose)
```
Process the volume to extract the ILM, RPE, ISOS and THICKNESS, save csv files in the destination folder with images of the scans and segmentation overimposed.
