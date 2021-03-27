#!/usr/bin/python3
import os
import sys
import argparse
import ROOT
import random

def main():
    SFLASH_FLUKA=os.environ.get('SFLASH_FLUKA')
    if SFLASH_FLUKA == None:
        sys.stderr.write('WARNING: SFLASH_FLUKA env variable not set\n')
        SFLASH_FLUKA=os.path.dirname(os.path.realpath(__file__))
        sys.stderr.write('Using %s\n' % SFLASH_FLUKA)
    sflashAnalysis_C=SFLASH_FLUKA+'/'+'sflashAnalysis.C'
    if not os.path.isfile(sflashAnalysis_C):
        sys.stderr.write('Error: %s not found\n' % (sflashAnalysis_C))
        exit(2)
    sflashAnalysis_C=os.path.abspath(sflashAnalysis_C)

    ROOT.gROOT.LoadMacro(sflashAnalysis_C+'+')

    parser = argparse.ArgumentParser()
    parser.add_argument('files',nargs='*',
                        help='pass .root file names w/o prefixes or switches')
    parser.add_argument('-i', action='store', dest='listfile',
                        help='give an ascii list file with paths to a bunch of .root files')
    parser.add_argument('-f', action='store_true',
                        default=False,
                        dest='fOverwriteMode',
                        help='Overwrite output files')
    parser.add_argument('-o', action='store', 
                        dest='outdir', default="./",
                        help='Output directory (output files generated automatically)')

    if (len(sys.argv)==1):
        sys.stdout.write("\n");
        sys.stdout.write("Run sflashAnalysis.C on a bunch of waveform .root files\n");
        parser.print_help()
        sys.stdout.write("\n\n")
        sys.exit(1)

    args = parser.parse_args()

    flist_rel=[]
    if args.files != None:
        flist_rel = args.files
    
    if args.listfile != None:
        with open(args.listfile,"r") as f:
            flist_rel.extend(map(lambda s: s.strip(), f.readlines()))
    
    flist=map(lambda s: os.path.abspath(s), flist_rel)

    if(len(flist) < 1):
        sys.stderr.write("No input files; stop\n")
        exit(0)
    if(not os.path.isdir(args.outdir)):
        sys.stderr.write("error: -o: output directory {:s} not found\n".format(args.outdir))
        exit(0)
    outdir=os.path.abspath(args.outdir)
    
    for infile in flist:
        b=os.path.basename(infile)
        b = b[:-len(".root")]+".sflashAnalysis.root" if b.endswith(".root") else b+".sflashAnalysis.root"
        outfile=outdir+"/"+b
        if os.path.isfile(outfile) and not args.fOverwriteMode:
            sys.stderr.write("error: {:s} exists; use -f to overwrite output files\n".format(outfile));
            exit(2)
        sys.stdout.write("{:s} {:s}\n".format(infile,outfile))
        ROOT.sflashAnalysis(infile,outfile)


if __name__== '__main__':
    main()
