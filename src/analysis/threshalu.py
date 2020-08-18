alus = set()

f = open("alu_list.txt", "r")
for line in f.readlines():
	cname = line.split(" ")[0]
	alus.add(cname)
f.close()

thresholds = set()

g = open("results.thresh", "r")
for line in g.readlines():
	arr = line.split("\t")
	if arr[0] in alus and arr[1][3:] == "43g":
		thresholds.add(float(arr[2]))
		print(line[:-1])
g.close()

print(thresholds)
