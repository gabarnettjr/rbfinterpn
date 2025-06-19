#!/mask/software/RHEL7.7/bin/python3

# Given full paths to two sets of points (A and B), this script writes an output
# file that contains the index of the nearest neighbors in set B to each point
# in set A.  The output file has the same number of rows as the set A file.

# Greg Barnett
# May 2024

import sys
import numpy as np
from scipy.spatial import KDTree

if len(sys.argv) != 5 :
    raise ValueError("Exactly four inputs (aFile, bFile, outFile, numNeighbors) are required.")

################################################################################

# Load the points from set A.

pathToPointsA = sys.argv[1]
count = 0

with open(pathToPointsA) as f :
    for line in f :
        count += 1

pointsA = np.zeros((count, 2), np.float64)
count = 0

with open(pathToPointsA) as f :
    for line in f :
        tmp = line.split(',')
        pointsA[count, 0] = np.float64(tmp[0].strip())
        pointsA[count, 1] = np.float64(tmp[1].strip())
        count += 1

# for i in range(count) :
    # print(str(pointsA[i, 0]) + ", " + str(pointsA[i, 1]))
# print()

################################################################################

# Load the points from set B.

pathToPointsB = sys.argv[2]
count = 0

with open(pathToPointsB) as f :
    for line in f :
        count += 1

pointsB = np.zeros((count, 2), np.float64)
count = 0

with open(pathToPointsB) as f :
    for line in f :
        tmp = line.split(',')
        pointsB[count, 0] = np.float64(tmp[0].strip())
        pointsB[count, 1] = np.float64(tmp[1].strip())
        count += 1

# for i in range(count) :
    # print(str(pointsB[i, 0]) + ", " + str(pointsB[i, 1]))
# print()

################################################################################

# Make the KDTree object and print nearest neighbor indices to the output file.

pathToOutput = sys.argv[3]
tree = KDTree(pointsB)
dst, ind = tree.query(pointsA, k = int(sys.argv[4]))

with open(pathToOutput, 'w') as f :
    for i in range(np.shape(ind)[0]) :
        for j in range(np.shape(ind)[1] - 1) :
            f.write("{0:10d}, ".format(ind[i, j]))
        f.write("{0:10d}\n".format(ind[i, -1]))
