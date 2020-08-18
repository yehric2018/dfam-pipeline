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

from sequence_util import nearestDivergence

def getAlus():
	f = open("alu_list.txt", "r")
	alus = []
	for line in f.readlines():
		alu = line.split()
		consensus = alu[0]
		div = nearestDivergence(float(alu[1].split("=")[1]))
		alus.append((consensus, str(div)))
	return alus

def convertToBed(fname, bed, threshold):
	print("converting")
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

def alusToBed(threshold, alus, bpath, dirpath):
	bed = open(bpath, "w")
	for alu in alus:
		fpath = dirpath + alu[0] + "/" + alu[0] + "_" + alu[1] + "p43g.sc"
		convertToBed(fpath, bed, threshold)
	bed.close()

def makeBeds(threshold, alus):
	alusToBed(threshold, alus, "genomic_alus.bed", "../../results/genomic_hits/")
	alusToBed(threshold, alus, "benchmark_alus.bed", "../../results/benchmark_hits/")

def splitChromosomes(fpath):
	f = open(fpath)
	for line in f.readlines():
		print()
	f.close()

def main():
	parser = argparse.ArgumentParser()
	args = parser.parse_args()

	alus = getAlus()

	makeBeds(120.05, alus)


if __name__ == '__main__':
    main()
