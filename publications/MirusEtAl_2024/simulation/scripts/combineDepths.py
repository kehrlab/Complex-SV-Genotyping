# Tim Mirus
# small python helper script to calculate the average coverage across two chromosomes
# given the depths and lengths for the two chromosomes

import sys

args = sys.argv

if len(args) != 5:
    print("Usage: combineDepths.py <depth1> <length1> <depth2> <length2>")
    exit(1)

d1 = float(args[1])
d2 = float(args[3])
l1 = float(args[2])
l2 = float(args[4])
print((d1*l1+d2*l2)/(l1+l2))
