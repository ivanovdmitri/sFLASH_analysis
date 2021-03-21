# sFLASH_analysis
Analysis of sFLASH data. sFLASH (Super FLASH, FLASH = FLuorescence in Air from SHowers) experiminent was carried out at the Stanford Linear Accelerator End Station A facility with an aim of measuring the fluorescence yield of electromangetic showers produced by a 14.49 GeV electron beam hitting Al2O3 targets.   

## conv_wf2rt.py
A python script that
* Reads all firmware output, parses csv files and combines all relevant records
  into events that are readouts of all sensors corresponding to each accelerator
  beam pulse
* Reads supplementary setting and calibration xlsx or csv files and attaches
  relevant information to each event
* Records events in ROOT tree format https://root.cern.ch/doc/master/group__tree.html.

## sflashAnalysis.C
C++ code that analyses sFLASH data prepared by runing conv_wf2rt.py on the raw data and settings files. The code 
* analyses the signal waveforms, recorded in nano-second time slices, calculates and subtracts the pedestals and determines the pulse areas 
* applies the calibration to the photomultiplier tube pulse areas and to the charge measuring coil
* calcualtes variables that are needed for applying event selection quality cuts
* contains routines that calculate coil-calibrated signal and noise levels, calculate the fluorescence yield for each target configuration, and plot the results 
