f = open("seq.txt", "r")
line = f.readline().lower()

cg = 0
tot = 0
for i in line:
    if i in "gcat":
        if i in "gc":
            cg += 1
        tot += 1

print(100.0 * cg / tot)
