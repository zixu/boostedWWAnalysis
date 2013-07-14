#! /usr/bin/env python
import os
import glob
import math
import array
import re

import subprocess
from subprocess import Popen

import sys

from optparse import OptionParser
parser = OptionParser()

parser.add_option('-i', '--inputfile',action="store",type="string",dest="inputfile",default="")
parser.add_option('-m', '--lumibefore',action="store",type="float",dest="lumibefore",default=1)
parser.add_option('-n', '--lumiafter' ,action="store",type="float",dest="lumiafter" ,default=1)


(options, args) = parser.parse_args()

if __name__ == '__main__':

    f=open(options.inputfile);
    file=f.readlines();

    lumiscale=options.lumiafter/options.lumibefore;

    for i in range(0,34):
        if re.match("observation",file[i]):
            line=file[i][11:-1];
            line=line.strip("\n");
            result=[];#
            for number in line.split(" "):
                if number: result.append(float(number)*lumiscale );
            print "observation %s"%(result[0])
        elif re.match("rate",file[i]):
            line=file[i][4:-1];
            result=[];#
            for number in line.split(" "):
                if number: result.append(float(number)*lumiscale );
            print "rate %s %s %s %s %s"%(result[0],result[1],result[2],result[3],result[4] )
        else:
            print file[i].strip("\n");

            

