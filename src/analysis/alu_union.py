#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
alu_union.py - Union all alu hits together for a given threshold and
find the TP/FP ratio.

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import argparse
import os
import re

def convertToBed(fname, bed, threshold):
	f = open(fname, "r")
	scoreRegex = re.compile(r"^\s*(\d+)(\s+\d+\.\d+){3}\s+(\d+):(\d+)\-\d+\s+(\d+)\s+(\d+)")

	line = f.readline()
	while line != "":
		mo = scoreRegex.search(line)
		if mo and int(mo.group(1)) > threshold:
			start = int(mo.group(4)) + int(mo.group(5))
			end = int(mo.group(4)) + int(mo.group(6))
			bed.write(mo.group(3) + "\t" + str(start) + "\t" + str(end) + "\n")
		line = f.readline()

def main():
	parser = argparse.ArgumentParser()
	args = parser.parse_args()

	threshold = 120.05
	bed = open("14p35g.bed", "w")
	
	convertToBed("../../results/genomic_hits/DF0000002/DF0000002_14p35g.sc", bed, threshold)
	convertToBed("../../results/genomic_hits/DF0000003/DF0000003_14p35g.sc", bed, threshold)

if __name__ == '__main__':
    main()
