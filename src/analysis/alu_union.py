import re
import os

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

def main():
	alus = getAlus()
	genomic_hits = getAllHits(alus, "../../results/genomic_hits/")
	print(genomic_hits)
	benchmark_hits = getAllHits(alus, "../../results/benchmark_hits/")
	print(benchmark_hits)

if __name__ == '__main__':
    main()
