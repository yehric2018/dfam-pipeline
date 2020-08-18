import re
import os
import subprocess

from sequence_util import nearestDivergence

def getAlus():
	f = open("test_list.txt", "r")
	alus = []
	for line in f.readlines():
		alu = line.split()
		consensus = alu[0]
		div = nearestDivergence(float(alu[1].split("=")[1]))
		alus.append((consensus, str(div)))
	return alus

def getHits(fname):
	hits = {}

	f = open(fname, "r")
	regex = re.compile(r"^\s*(\d+)(\s+\d+\.\d+){3}\s+(\d+):(\d+)\-\d+\s+(\d+)\s+(\d+)")

	line = f.readline()
	while line != "":
		mo = regex.search(line)
		if mo:
			chrom = mo.group(3)
			start = int(mo.group(4)) + int(mo.group(5))
			end = int(mo.group(4)) + int(mo.group(6))
			score = int(mo.group(1))
			if chrom not in hits:
				hits[chrom] = []
			hits[chrom].append(((start, end), score))
		line = f.readline()
	f.close()

	return hits

def getAllHits(alus, path):
	hits = {}
	for alu in alus:
		fname = alu[0] + "/" + alu[0] + "_" + alu[1] + "p43g.sc"
		new_hits = getHits(os.path.join(path, fname))
		for key in new_hits:
			if key not in hits:
				hits[key] = []
			hits[key] = hits[key] + new_hits[key]
	return hits

def sort1(elem):
	return elem[0][1]

def sort2(elem):
	return elem[0][0]

def writeToBed(chrom, hits_list, threshold):
	# Filter hits
	hits = list(filter(lambda x : x[1] > threshold, hits_list))
	print(len(hits))

	# Takes in a list of hits on one chromosome
	# Sorts the list in this function
	hits.sort(key=sort1)
	hits.sort(key=sort2)
	bed = open(chrom + str(threshold) + ".bed", "w") # temp.bed
	for hit in hits:
		bed.write(chrom + "\t" + str(hit[0][0]) + "\t" + str(hit[0][1]) + "\n")
	bed.close()

def countMerged(hits, threshold):
	# hits is a dict, chromosomes are keys and list of hits are vals
	count = 0
	for chrom in hits:
		print("running on chromosome " + chrom + ", threshold " + str(threshold))
		print(len(hits[chrom]))
		writeToBed(chrom, hits[chrom], threshold)
		ps = subprocess.Popen((["bedtools", "merge", "-i", "temp.bed"]), stdout=subprocess.PIPE)
		output = subprocess.check_output(("wc", "-l"), stdin=ps.stdout)
		ps.wait()

		count += int(output.strip())
		print(int(output.strip()))
	return count

def getFDR(genomic_hits, benchmark_hits, threshold):
	gcount = countMerged(genomic_hits, threshold)
	bcount = countMerged(benchmark_hits, threshold)
	fdr = 1.0 * bcount / gcount
	print(str(threshold) + "\t" + str(gcount) + "\t" + str(bcount) + "\t" + str(fdr))
	return fdr

def main():
	alus = getAlus()
	genomic_hits = getAllHits(alus, "../../results/genomic_hits/")
	benchmark_hits = getAllHits(alus, "../../results/benchmark_hits/")

	thresholds = [128.05, 129.05, 130.05, 131.05, 134.05, 123.05, 125.05, 126.05, 127.05]
	for thresh in thresholds:
		fdr = getFDR(genomic_hits, benchmark_hits, thresh)

if __name__ == '__main__':
    main()
