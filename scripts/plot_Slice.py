#!/usr/bin/env python3
# coding=UTF-8

import sys, time
import os.path
import argparse
import numpy as np
from TimeDomainField import *
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the input file')
parser.add_argument('-ix', help="index x for slicing", dest="ix", type=int)
parser.add_argument('-iy', help="index y for slicing (default -ix=Nx/2)", dest="iy", type=int)
parser.add_argument('-Emax', help="plot range for the E field [V/m]", dest="emax", type=float)
parser.add_argument('-Bmax', help="plot range for the B field [T]", dest="bmax", type=float)

print()
args = parser.parse_args()

infile = args.file
fileOK = os.path.isfile(infile)
if not fileOK:
  print("file not found")
  sys.exit()

Source = TimeDomainField()
# Open the file for reading
print("reading %s" % infile)
Source.read(infile)

print("Nx=%d Ny=%d" % (Source.Nx,Source.Ny))
print("t0=%g dt=%g NOTS=%d" % (Source.t0, Source.dt, Source.Nt))

if args.iy != None:
    iy = args.iy
    print("slicing iy=%d",iy)
else:
    if args.ix != None:
        ix = args.ix
    else:
        ix = Source.Nx // 2
    print("slicing ix=%d",ix)

