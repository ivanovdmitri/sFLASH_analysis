# sFLASH_analysis
Analysis of sFLASH experiminent carried out at Stanford Linear Accelerator in
2018.

## conv_wf2rt.py
A python script that
* Reads all firmware output, parses csv files and combines all relevant records
  into events that are readouts of all sensors corresponding to each accelerator
  beam pulse
* Reads supplementary setting and calibration xlsx or csv files and attaches
  relevant information to each event
* Records events in ROOT tree format https://root.cern.ch/doc/master/group__tree.html.
