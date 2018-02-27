#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################################
#     calculate phasiRNA abundance       #
#      from featureCount result          #
#             2017.12.12                 #
##########################################
__author__ = "K.R.Chow"

import sys
import os
import argparse
import re
from multiBioPro import MultiSys, FeatureCount
from collections import defaultdict


def Version(options, script, version):
    if options.version:
        print("{0}: {1} version".format(script, version))


script = os.path.basename(sys.argv[0])
version = 'v1.0'

parser = argparse.ArgumentParser(
    description='This is a Python3-based program for calculatting phasiRNA abundance from featureCount result.')
parser.add_argument('-i', '--input', action="store",
                    help='The featureCount counts')
parser.add_argument('-m', '--mode', action="store",
                    default="rpm", help='The featureCount counts')
parser.add_argument('-n', '--name', action="store", nargs='+',
                    default="rpm", help='The featureCount counts')
parser.add_argument('-o', '--output', action="store",
                    default="result.txt", help='The output result')
parser.add_argument('-v', '--version', action='store_true',
                    default=False, help='The srcipt version')
MultiSys.ArgparseHelp(parser)
options = parser.parse_args()
Version(options, script, version)

inputFile = options.input
mode = options.mode
nameList = options.name
outputFile = options.output
MultiSys.FileExist(inputFile)

lineList = FeatureCount.Exp(file=inputFile, mode=mode)
if nameList:
    lineList[0][6:] = nameList

with open(outputFile, 'w') as f:
    f.write('#Expression level:'+mode + '\n')
    for line in lineList:
        f.write('\t'.join(MultiSys.List2Str(line)) + '\n')
