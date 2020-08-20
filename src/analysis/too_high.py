import re

# For score thresholds that are too high!
# Find the low complexity sequences that score high in benchmark genome

def main(fpath):
	consensus_size = 5102
	counts = [0] * consensus_size
	ranges = []

	regex1 = re.compile(r"(\d+)\s+(\d+)\s+\(\d+\)\s*$")
	regex2 = re.compile(r"\(\d+\)\s+(\d+)\s+(\d+)\s*$")
	
	f = open("../../results/benchmark_hits/DF0000187/DF0000187_" + fpath, "r")
	line = f.readline()
	while line != "":
		mo1 = regex1.search(line)
		mo2 = regex2.search(line)
		if mo1:
			start = int(mo1.group(1))
			end = int(mo1.group(2))
			ranges.append((start, end))
		if mo2:
			start = int(mo2.group(2))
			end = int(mo2.group(1))
			ranges.append((start, end))
		line = f.readline()
	f.close()

	for r in ranges:
		for x in range(r[0], r[1]):
			counts[x] += 1

	return counts
	
if __name__ == '__main__':
	seqs = []
	seqs.append(main("14p35g.sc"))
	seqs.append(main("14p37g.sc"))
	seqs.append(main("14p39g.sc"))
	seqs.append(main("14p41g.sc"))
	seqs.append(main("14p43g.sc"))
	seqs.append(main("14p45g.sc"))
	seqs.append(main("14p47g.sc"))
	seqs.append(main("14p49g.sc"))
	seqs.append(main("14p51g.sc"))
	seqs.append(main("14p53g.sc"))

	final = [0] * 5102
	for x in range(5102):
		final[x] = seqs[0][x] + seqs[1][x] + seqs[2][x] + seqs[3][x] + seqs[4][x] + seqs[5][x] + seqs[6][x] + seqs[7][x] + seqs[8][x] + seqs[9][x]
	print(final)
