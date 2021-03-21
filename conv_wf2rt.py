#!/usr/bin/env python

import argparse
from ROOT import TFile, TTree
from array import array
import itertools
import numpy as np
import sys
import os
import re


def read_xlsx_sheet(xlsx_file,sheet_name):
    # one can directly read excel files if
    # python-xlrd package has been installed 
    try:
        import xlrd
    except ImportError:
        sys.stderr.write("ERROR: package python-xlrd not intalled ")
        sys.stderr.write("but an .xlsx file \n\'{:s}\' given!\n".\
                         format(xlsx_file))
        sys.stderr.write("Either convert it to .csv format or ")
        sys.stderr.write("install python-xlrd on your system!\n")
        sys.exit(2)
    with xlrd.open_workbook(xlsx_file) as wb:
        try:
            sheet=wb.sheet_by_name(sheet_name)
        except:
            sys.stderr.write("ERROR: no sheet named \'{:s}\' in file {:s}!\n".\
                             format(sheet_name,xlsx_file))
            sys.exit(2)
        if sheet.nrows < 1:
            sys.stderr.write("ERROR: sheet \'{:s}\' in file {:s} has no rows!\n".\
                             format(sheet_name,xlsx_file))
            sys.exit(2)
        header=map(lambda s: str(s.encode("ascii")), sheet.row_values(0))
        values=[]
        for row in range(1,sheet.nrows):
            row_values={}
            for col in range(0,sheet.ncols):
                v=sheet.cell(row,col).value
                if type(v) == unicode:
                    v = str(v.encode("ascii"))
                row_values[header[col]] = v
            values.append(row_values)
    return (header, values)
            
def read_csv_file(csv_file):
    # one can directly read CSV files if
    # python-csv package has been installed 
    try:
        import csv
    except ImportError:
        sys.stderr.write("ERROR: package python-csv not intalled ")
        sys.stderr.write("but a .csv file {:s} given! \nFigure out ".format(csv_file))
        sys.stderr.write("how to install python-csv package on your system !\n")
        sys.exit(2)
    header=[]
    values=[]
    with open(csv_file) as csvfile:
        reader = csv.reader(csvfile)
        header_read = next(reader,None)
        header.extend(header_read)
    with open(csv_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            row_values={}
            for name,val in row.iteritems():
                if type(val) == unicode:
                    val = str(val.encode("ascii"))
                row_values[name]=val
            values.append(row_values)
    return (header,values)


class settings_class:
    def __init__(self):
        # Variables that should be present in run settings CVS file
        # (e.g. sFLASH_run3_settings.cvs that's made from
        # sFLASH_run3_settings.xlsx)
        self.SCOPE=int(0)          # oscilloscope ID
        self.CHANNEL=int(0)        # channel ID of the oscilloscope
        self.PMT_SERIAL=str("")    # PMT serial code, if applies, or COIL
                                   # if the coil is what's plugged in to the
                                   # scope channel
        self.CONNECTED_TO=str("")  # brief description of what's plugged in
                                   # to the scope channel 
        self.REMARKS=str("")       # additional remarks about what's plugged
                                   # in to the scope channel

def load_settings(settings_file):
    sys.stdout.write("loading settings file {:s} ...\n".format(settings_file))
    if settings_file.endswith(".csv"):
        (header,values) = read_csv_file(settings_file)
    elif settings_file.endswith(".xlsx"):
        (header,values) = read_xlsx_sheet(settings_file,"settings")
    else:
        sys.stderr.write("ERROR: file {:s} doesn\'t end with .csv or .xlsx !".format(settings_file));
        sys.exit(2)
    settings = []
    for row in values:
        pmti=settings_class()
        for name,val in row.iteritems():
            vars(pmti)[name] = type(vars(pmti)[name])(val)
        settings.append(pmti)
    sys.stdout.write("settings file {:s} loaded successfully\n".format(settings_file))
    sys.stdout.flush()
    return (header,settings)

class conditions_class:
    def __init__(self):
        # Variables that should be present in the CSV conditions file
        # (E.G. sFLASH_run3_conditions.csv is made from sFLASH_run3_conditions.xlsx by using SaveAs (CVS format) of MS Excel or OpenOffice Spreadsheet)
        self.yyyymmdd_start=int(0) # date+time at the sub-run start
        self.hhmmss_start=int(0)
        self.yyyymmdd_end=int(0)   # date+time at the sub-run end
        self.hhmmss_end=int(0)
        self.run_id=int(0)         # ID of the sub-run
        self.run_type=int(0)       # type of the sub-run (beam or UVLED ?)
        self.rl=int(0)             # number of radiation lenghts
        self.nevent=int(0)         # number of events
        self.blind=int(0)          # blind open or closed?
        self.shutter=int(0)        # shutter open or closed?
        self.status=int(0)         # sub-run status, good or bad?
        self.E_GeV=float(0)        # Beam energy in GeV
        self.C_per_Vs=float(0)     # COIL waveform integral to charge converter in Coulombs per (Voltage * Second)
        self.PMT_1_HV=float(0)     # PMT high voltage
        self.PMT_2_HV=float(0)
        self.PMT_3_HV=float(0)
        self.PMT_4_HV=float(0)
        self.PMT_A_HV=float(0)
        self.PMT_B_HV=float(0)
        self.PMT_C_HV=float(0)

        # Parameters to be used in PMT calibration formula: 
        # Gain (in NPE_per_Vs) = 1.0 / [ EXP(PMT_X_ALPHA)*(PMT_X_HV-PMT_X_OFFSET)^PMT_X_BETA
        self.PMT_1_OFFSET=float(0) # OFFSET in the PMT calibration formula
        self.PMT_1_ALPHA=float(0)  # ALPHA parameter for the PMT calibration formula
        self.PMT_1_BETA=float(0)   # BETA parameter for the PMT calibration formula
        self.PMT_2_OFFSET=float(0)
        self.PMT_2_ALPHA=float(0)
        self.PMT_2_BETA=float(0)
        self.PMT_3_OFFSET=float(0)
        self.PMT_3_ALPHA=float(0)
        self.PMT_3_BETA=float(0)
        self.PMT_4_OFFSET=float(0)
        self.PMT_4_ALPHA=float(0)
        self.PMT_4_BETA=float(0)
        self.PMT_A_OFFSET=float(0)
        self.PMT_A_ALPHA=float(0)
        self.PMT_A_BETA=float(0)
        self.PMT_B_OFFSET=float(0)
        self.PMT_B_ALPHA=float(0)
        self.PMT_B_BETA=float(0)
        self.PMT_C_OFFSET=float(0)
        self.PMT_C_ALPHA=float(0)
        self.PMT_C_BETA=float(0)
        
        # Energy tracing factors for active regions up to the blind (TB) and up to the shutter (TS)
        # as calculated by FLUKA and GEANT4 simulations and the corresponding raytracing routines.
        self.PMT_1_TB_FK=float(0)  # Energy tracing factor for the active region until blind calculated by FLUKA simulation, [MeV / pC]
        self.PMT_2_TB_FK=float(0)  
        self.PMT_3_TB_FK=float(0)
        self.PMT_4_TB_FK=float(0)
        self.PMT_A_TB_FK=float(0)
        self.PMT_B_TB_FK=float(0)
        self.PMT_C_TB_FK=float(0)
        self.PMT_1_TS_FK=float(0)  # Energy tracing factor for the active region until shutter calculated by FLUKA simulation, [MeV / pC]
        self.PMT_2_TS_FK=float(0)  
        self.PMT_3_TS_FK=float(0)
        self.PMT_4_TS_FK=float(0)
        self.PMT_A_TS_FK=float(0)
        self.PMT_B_TS_FK=float(0)
        self.PMT_C_TS_FK=float(0)
        self.PMT_1_TB_G4=float(0)  # Energy tracing factor for the active region until blind calculated by Geant4 simulation, [MeV / pC]
        self.PMT_2_TB_G4=float(0)  
        self.PMT_3_TB_G4=float(0)
        self.PMT_4_TB_G4=float(0)
        self.PMT_A_TB_G4=float(0)
        self.PMT_B_TB_G4=float(0)
        self.PMT_C_TB_G4=float(0)
        self.PMT_1_TS_G4=float(0)  # Energy tracing factor for the active region until shutter calculated by Geant4 simulation, [MeV / pC]
        self.PMT_2_TS_G4=float(0)  
        self.PMT_3_TS_G4=float(0)
        self.PMT_4_TS_G4=float(0)
        self.PMT_A_TS_G4=float(0)
        self.PMT_B_TS_G4=float(0)
        self.PMT_C_TS_G4=float(0)


def load_conditions(conditions_file):
    sys.stdout.write("loading conditions file {:s} ...\n".format(conditions_file))
    if conditions_file.endswith(".csv"):
        (header,values) = read_csv_file(conditions_file)
    elif conditions_file.endswith(".xlsx"):
        (header,values) = read_xlsx_sheet(conditions_file,"conditions")
    else:
        sys.stderr.write("ERROR: file {:s} doesn\'t end with .csv or .xlsx !".format(conditions_file));
        sys.exit(2)
    conditions = {}
    for row in values:
        cond=conditions_class()
        for name,val in row.iteritems():
            vars(cond)[name] = type(vars(cond)[name])(val)
        conditions[cond.run_id] = cond
    sys.stdout.write("conditions file {:s} loaded successfully\n".format(conditions_file))
    sys.stdout.flush()
    return (header,conditions)



class wf_class:
    def __init__(self):
        # class that describes the waveform data in RunNNNNNN-YYYYMMDD-HHMMSS-wf-N.csv type of file
        self.t0=float(0)
        self.seq=int(0)
        self.dt=[float(0),float(0),float(0),float(0)]
        self.data=[[],[],[],[]]
    def __str__(self):
        return '{0:.3f} {1:d} {2:.6e} {3:.6e} {4:.6e} {5:.6e} {6:d}'\
            .format(self.t0,self.seq,self.dt[0],\
                    self.dt[1],self.dt[2],self.dt[3],len(self.data[0]))

def parse_wf_file_name(fname):
    pattern_run_header=re.compile(r"""
.*Run(?P<run_id>\d{6})
-(?P<yyyymmdd>\d{8})-(?P<hhmmss>\d{6})-wf-(?P<scope_id>\d).csv
""",re.VERBOSE)
    match=pattern_run_header.match(fname)
    if match == None:
        raise Exception("{:s} is not a RunNNNNNN-YYYYMMDD-HHMMSS-wf-N.csv type file".format(fname))
    run_id = int(match.group("run_id"))
    yyyymmdd=int(match.group("yyyymmdd"))
    hhmmss=int(match.group("hhmmss"))
    scope_id=int(match.group("scope_id"))
    return (run_id,yyyymmdd,hhmmss,scope_id)

def parse_waveform_data_file(infile): 
    sys.stdout.write("parsing waveform data file {:s} ...\n".format(infile))
    Waveforms=[]
    lines=map(lambda s: s.strip(), open(infile,"r").readlines())
    wf=wf_class()
    i=0
    while i < len(lines):
        line=lines[i]
        if line == "t0,seq":
            
            # header of a new waveform.  store the previously filled waveform if it
            # has been filled
            if(len(wf.data[0])>0):
                Waveforms.append(wf)
            
            # check if got an empty header at the end of the file
            if i+5 > len(lines):
                sys.stderr.write("warning: empty waveform (no waveform data) at the end of the file:\n");
                for k in range(i,len(lines)):
                    sys.stderr.write("{:s} line {:d}: {:s}\n".format(infile,k,lines[k]))
                i= len(lines)
                continue
            
            # store header lines
            header_lines=lines[i:i+5]
            
            # check if this is an emtpy waveform or not
            wf_empty=False
            for j in range(1,5):
                if (header_lines[j] == "t0,seq"):
                    sys.stderr.write("warning: empty waveform (no waveform data):\n");
                    for k in range(0,j+1):
                        sys.stderr.write("{:s} line {:d}: {:s}\n".format(infile,i+k,lines[i+k]))
                    i=i+j
                    wf_empty=True
                    break
            if wf_empty:
                continue
                
            # make sure that the header isn't corrupted.  if it is corrupted,
            # continue reading until a new header is found or the end of the file
            # is reached.
            corrupted_header=False
            if header_lines[4] != "WaveForm Data":
                sys.stderr.write("warning: corrupted header:\n");
                for k in range(0,5):
                    sys.stderr.write("{:s} line {:d}: {:s}\n".format(infile,i+k,lines[i+k]))
                i=i+5
                corrupted_header=True
                break
            if corrupted_header:
                while (i < len(lines) and  lines[i] != "t0,seq"):
                    i=i+1
                continue # be sure to re-iterate all the above code on the new header
            
            # start a new waveform
            wf=wf_class()
            
            # fill the header information for the new waveform
            wf.t0=float(header_lines[1].split(',')[0])
            wf.seq=int(header_lines[1].split(',')[1])
            for k in range(0,4):
                wf.dt[k]=float(header_lines[3].split(',')[k])
            wf.data=[[],[],[],[]]
            i=i+5
            continue
        l=line.split(',')
        for k in range(0,4):
            wf.data[k].append(float(l[k]))
        i=i+1
        if(i == len(lines)):
            Waveforms.append(wf)
            wf=None
    sys.stdout.write("done parsing waveform data file {:s} \n".format(infile))
    sys.stdout.flush();
    return Waveforms


def convert_waveform_files(header_and_settings,wf0_infile,wf1_infile,header_and_conditions,outfile):

    def make_rt_branch(name,vals_read,n_item_variable_name=""):
        # take the name of the variables, the list of values,
        # and a variable (stored in ROOT tree) that describes
        # the number of item, (if empty, assume the number of entries is 1)
        # and prepare all necessary components of a ROOT - tree branch:
        # (branch name, branch container, branch descriptor)
        # and return this as a list of the three things above
        if(len(vals_read) > 1 and len(n_item_variable_name) == 0):
            sys.stderr.write("ERROR: make_rt_branch: number of values passed greater than 1 but no")
            sys.stderr.write("root tree variable name given that describes the number of items!\n")
            sys.exit(2)

        n_item_descriptor="";
        if len(vals_read) > 1:
            n_item_descriptor = "["+n_item_variable_name+"]";

        if type(vals_read[0]) == str:
            vals_adj=map(lambda s : s+"\0", vals_read)
            l_max=max(map(lambda s: len(s), vals_adj))
            vals=map(lambda s : s+(l_max-len(s))*"\0", vals_adj)
            vals_array=np.arange(len(vals)*l_max,dtype='b').reshape((len(vals),l_max),order="C")
            i=0
            for val in vals:
                j=0
                for v in val:
                    vals_array[i][j] = ord(v)
                    j=j+1
                i=i+1
            return (name,vals_array,name+n_item_descriptor+"["+str(l_max)+"]/B")
        elif type(vals_read[0]) == int:
            vals_array=array("i")
            vals_array.extend(vals_read)
            return (name,vals_array,name+n_item_descriptor+"/I")
        elif type(vals_read[0]) == float:
            vals_array=array("d")
            vals_array.extend(vals_read)
            return (name,vals_array,name+n_item_descriptor+"/D")
        else:
            sys.stderr.write("ERROR: unsupported data type!\n")
            sys.exit(2)

    # settings variables that apply to all run_id's of the global run
    settings_rt_branches=[]
    if len(header_and_settings) == 2:
        (header,settings) = header_and_settings
        if(len(settings) > 0):
            # ROOT tree variable that describes the number of condition information entries
            settings_rt_branches.append(make_rt_branch("n_connected",[int(len(settings))]))
            for name in header:
                vals_read=map(lambda s: vars(s)[name], settings)
                settings_rt_branches.append(make_rt_branch(name,vals_read,settings_rt_branches[0][0]))


    # waveforms from both oscilloscopes for the particular sub-run,
    # labeled by the run_id variable. Must have two waveform files,
    # from the two oscilloscopes, that have the same run_id
    waveform_static_rt_branches=[]
    Waveforms0 = parse_waveform_data_file(wf0_infile)
    Waveforms1 = parse_waveform_data_file(wf1_infile)
    nwf0=len(Waveforms0)
    nwf1=len(Waveforms1)
    nwf=nwf0
    if(nwf0 != nwf1):
        sys.stderr.write("WARNING: number of waveforms in file {:s} {%d} is not the same\n".format(wf0_infile,nwf0));
        sys.stderr.write("as the number of waveforms in file {:s} {:d}!\n".format(wf1_infile,nwf1));
        if nwf1 < nw0:
            nwf = nwf1
    (run_id,yyyymmdd,hhmmss,scope_id) = parse_wf_file_name(wf0_infile)
    waveform_static_rt_branches.append(make_rt_branch("run_id",[int(run_id)]))
    waveform_static_rt_branches.append(make_rt_branch("yyyymmdd",[int(yyyymmdd)]))                                  
    waveform_static_rt_branches.append(make_rt_branch("hhmmss",[int(hhmmss)]))

    # look at the parsed conditions variables listed for each sub-run
    # pick out a set of conditions that's relevant for the waveforms being parsed, 
    # i.e. a set of conditions that correspond to the run_id of the waveforms
    conditions_rt_branches=[]
    if len(header_and_conditions) == 2:
        (header,conditions) = header_and_conditions
        if run_id in conditions.keys():
            relevant_conditions=[conditions[run_id]]
            for name in header:
                vals_read=map(lambda s: vars(s)[name], relevant_conditions)
                conditions_rt_branches.append(make_rt_branch(name,vals_read))

                
    # Initialize ROOT file and allocate the ROOT tree linked to
    # that file
    f = TFile(outfile,"recreate")
    if f.IsZombie():
        exit(2)
    t = TTree("tsFLASHwf","Tree with SFLASH 2018 Waveform Data")


    # set the branches for the overall run settings, if the run settings data is available
    # as well the condition and simulation branches for the particular sub run, 
    # mapped out by the run_id of the waveforms, if the corresponding run conditions data is available
    for rt_branch in itertools.chain(settings_rt_branches,conditions_rt_branches,waveform_static_rt_branches):
        t.Branch(*rt_branch)
        
    # set the wavform branches for the particular sub run mapped out by the run_id
    # variables of these branches may vary from pulse to pulse
    ntmax=5000
    t00   = array("d",[0.0])
    seq0  = array("i",[0])
    t01   = array("d",[0.0])
    seq1  = array("i",[0])
    nt0   = array("i",[0])
    nt1   = array("i",[0])
    dt_coil = array("d",[0.0])
    dt_pmt_1 = array("d",[0.0])
    dt_pmt_4 = array("d",[0.0])
    dt_pmt_c = array("d",[0.0])
    dt_pmt_2 = array("d",[0.0])
    dt_pmt_3 = array("d",[0.0])
    dt_pmt_a = array("d",[0.0])
    dt_pmt_b = array("d",[0.0])
    wfti0   = array("i",ntmax*[0])
    wfti1   = array("i",ntmax*[0])
    wf_coil = array("d",ntmax*[0.0])
    wf_pmt_1 = array("d",ntmax*[0.0])
    wf_pmt_4 = array("d",ntmax*[0.0])
    wf_pmt_c = array("d",ntmax*[0.0])
    wf_pmt_2 = array("d",ntmax*[0.0])
    wf_pmt_3 = array("d",ntmax*[0.0])
    wf_pmt_a = array("d",ntmax*[0.0])
    wf_pmt_b = array("d",ntmax*[0.0])
 
       
    t.Branch("t00",t00,"t00/D")
    t.Branch("seq0",seq0,"seq0/I")
    t.Branch("t01",t01,"t01/D")
    t.Branch("seq1",seq1,"seq1/I")
    t.Branch("dt_coil",dt_coil,"dt_coil/D")
    t.Branch("dt_pmt_1",dt_pmt_1,"dt_pmt_1/D")
    t.Branch("dt_pmt_4",dt_pmt_4,"dt_pmt_4/D")
    t.Branch("dt_pmt_c",dt_pmt_c,"dt_pmt_c/D")
    t.Branch("dt_pmt_2",dt_pmt_2,"dt_pmt_2/D")
    t.Branch("dt_pmt_3",dt_pmt_3,"dt_pmt_3/D")
    t.Branch("dt_pmt_a",dt_pmt_a,"dt_pmt_a/D")
    t.Branch("dt_pmt_b",dt_pmt_b,"dt_pmt_b/D")
    t.Branch("nt0",nt0,"nt0/I")
    t.Branch("nt1",nt1,"nt1/I")
    t.Branch("wfti0",wfti0,"wfti0[nt0]/I")
    t.Branch("wf_coil",wf_coil,"wf_coil[nt0]/D")
    t.Branch("wf_pmt_1",wf_pmt_1,"wf_pmt_1[nt0]/D")
    t.Branch("wf_pmt_4",wf_pmt_4,"wf_pmt_4[nt0]/D")
    t.Branch("wf_pmt_c",wf_pmt_c,"wf_pmt_c[nt0]/D")
    t.Branch("wfti1",wfti1,"wfti1[nt1]/I")
    t.Branch("wf_pmt_2",wf_pmt_2,"wf_pmt_2[nt1]/D")
    t.Branch("wf_pmt_3",wf_pmt_3,"wf_pmt_3[nt1]/D")
    t.Branch("wf_pmt_a",wf_pmt_a,"wf_pmt_a[nt1]/D")
    t.Branch("wf_pmt_b",wf_pmt_b,"wf_pmt_b[nt1]/D")

 
        
    # combine the waveforms into event and fill the tree
    for iwf in range(0,nwf):
        wf0        = Waveforms0[iwf]
        t00[0]     = wf0.t0
        seq0[0]    = wf0.seq
        dt_coil[0] = wf0.dt[0]
        dt_pmt_1[0] = wf0.dt[1]
        dt_pmt_4[0] = wf0.dt[2]
        dt_pmt_c[0] = wf0.dt[3]
        nt0[0]     = len(wf0.data[0])
        for i in range(0,nt0[0]):
            wfti0[i]   = i
            wf_coil[i] = wf0.data[0][i]
            wf_pmt_1[i] = wf0.data[1][i]
            wf_pmt_4[i] = wf0.data[2][i]
            wf_pmt_c[i] = wf0.data[3][i]
        wf1        = Waveforms1[iwf]
        t01[0]     = wf1.t0
        seq1[0]    = wf1.seq
        dt_pmt_2[0] = wf1.dt[0]
        dt_pmt_3[0] = wf1.dt[1]
        dt_pmt_a[0] = wf1.dt[2]
        dt_pmt_b[0] = wf1.dt[3]
        nt1[0]     = len(wf1.data[0])
        for i in range(0,nt1[0]):
            wfti1[i]   = i
            wf_pmt_2[i] = wf1.data[0][i]
            wf_pmt_3[i] = wf1.data[1][i]
            wf_pmt_a[i] = wf1.data[2][i]
            wf_pmt_b[i] = wf1.data[3][i]
        t.Fill()
    t.Write()
    f.Close()


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("files",nargs="*",
                        help="pass RunNNNNNN-YYYYMMDD-HHMMSS-wf-N.csv file names w/o prefixes or switches,"+ 
                        "and the number of files must be even because there are two oscilloscopes")
    parser.add_argument("-i", action="store", dest="listfile",
                        help=' <string> give an ascii list file with paths to a bunch of RunNNNNNN-YYYYMMDD-HHMMSS-wf-N.csv files (even number of files)')
    parser.add_argument("--tty", action='store_true', 
                        default=False, 
                        dest="tty_input",
                        help="pipe RunNNNNNN-YYYYMMDD-HHMMSS-wf-N.csv file names from stdin (even number of files)")
    parser.add_argument("-settings", action='store',dest='settings_file',default="sFLASH_run3_settings.xlsx",
                        help="<string> Give a .csv or .xlsx file with scope/channel/PMT settings for the entire run period"+\
                        " (sFLASH_run3_settings.xlsx)")
    parser.add_argument('-conditions', action='store',dest='conditions_file',default="sFLASH_run3_conditions.xlsx",
                        help="<string> Give a .csv or .xlsx file with conditions for all sub-runs of run3 (sFLASH_run3_conditions.xlsx)")
    parser.add_argument("-o", action="store", dest="outdir", default="./",
                        help="<string >specify the output directory for the root tree files")
    parser.add_argument("-good", action="store_true", dest="good_runs_only", default=False,
                        help="Use this flag to parse only good runs (according to the conditions file) and skip all others\n")
    
    if (len(sys.argv)==1):
        sys.stdout.write("\n");
        sys.stdout.write("Convert sFLASH RUN3 (November 2018 run) waveform files, settings, and conditions, into root trees, last update October 2019\n")
        sys.stdout.write("DI <dmiivanov@gmail.com>\n")
        parser.print_help()
        sys.stdout.write("\n\n")
        sys.exit(2)
    args = parser.parse_args()
    infiles_rel=[]
    if args.files != None:
        infiles_rel.extend(args.files)
    if args.listfile != None:
        with open(args.listfile,"r") as f:
            infiles_rel.extend(map(lambda s: s.strip(), f.readlines()))
    if len(infiles_rel) < 1:
        sys.stderr.write("No input files\n")
        sys.exit(2)
    for infile in infiles_rel:
        if not os.path.isfile(infile):
            sys.stderr.write("ERROR: {0:s} file not found\n".format(infile))
            sys.exit(2)
    infiles=map(lambda s: os.path.abspath(s), infiles_rel)
    outdir=str(args.outdir).rstrip('/')
    if not os.path.isdir(outdir):
        sys.stdout.write("ERROR: output directory doesn\'t exist!\n");
        sys.exit(2)


    # run settings
    if os.path.isfile(args.settings_file):
        HeaderAndSettings = load_settings(os.path.abspath(args.settings_file))
    else:
        sys.stderr.write("WARNING: run3 settings file {:s} not found, parsing data\n".format(args.settings_file))
        sys.stderr.write("without knowing which PMT/COIL is connected to which scope!\n")
        HeaderAndSettings = []
    
    # prepare waveform infile pairs
    infile_pairs_all={}
    infile_pairs_good={}
    for infile in infiles:
        (run_id,yyyymmdd,hhmmss,scope_id) = parse_wf_file_name(infile)
        ind=run_id*100000000*1000000+yyyymmdd*1000000+hhmmss
        if ind not in infile_pairs_all.keys():
            infile_pairs_all[ind] = {}
        infile_pairs_all[ind][scope_id] = infile
    for ind,fpair in infile_pairs_all.iteritems():
        if len(fpair) < 2:
            sys.stderr.write("warning: run_id {:d} yyyymmdd {:d} hhmmss {:d}: number of waveform files less than 2!\n".
                             format(ind//100000000//1000000,(ind%(100000000*1000000))//1000000,ind%1000000))
            continue
        if len(fpair) > 2:
            sys.stderr.write("warning: run_id {:d} yyyymmdd {:d} hhmmss {:d}: number of waveform files greater than 2!\n".
                             format(ind//100000000//1000000,(ind%(100000000*1000000))//1000000,ind%1000000))
            continue
        infile_pairs_good[ind] = fpair

    # run conditions
    if os.path.isfile(args.conditions_file):
        HeaderAndConditions = load_conditions(os.path.abspath(args.conditions_file))
    else:
        if args.good_runs_only:
            sys.stderr.write("ERROR: requested parsing only good runs but conditions file\n")
            sys.stderr.write("that determines which runs are good is absent!\n")
            sys.exit(2)
        sys.stderr.write("WARNING: run3 conditions file {:s} not found, parsing data\n".format(args.conditions_file))
        sys.stderr.write("without any conditions and calibration information!\n")
        HeaderAndConditions=[]

        
    if(len(infile_pairs_good) < 1):
        sys.stderr.write("WARNING: don\'t have any good pairs of waveform files to analyze\n")
        

    # convert waveform data files to root trees
    for ind,fpair in infile_pairs_good.iteritems():
        wf0_infile=fpair[0]
        wf1_infile=fpair[1]
        # if parsing of only good runs is requested
        if args.good_runs_only:
            # skip the run if it's not in the list of runs
            # for which conditions are available (conditions should be
            # avaialble for all good runs)
            run_id=ind//100000000//1000000
            if not HeaderAndConditions[1].has_key(run_id):
                continue
            # skip the run if its status is not 1 (good run)
            if HeaderAndConditions[1][run_id].status != 1:
                continue
        outfile=outdir+"/"+os.path.basename(wf0_infile)
        outfile=outfile.replace("-0.csv",".root")
        convert_waveform_files(HeaderAndSettings,wf0_infile,wf1_infile,HeaderAndConditions,outfile)

    sys.stdout.write("\nDone\n")
        

if __name__ == "__main__":
    main()
