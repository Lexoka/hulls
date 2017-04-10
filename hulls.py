#!/usr/bin/env python3

import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

ANGLES				= [0, 30, 60, 90, 120, 180]
FREQUENCIES			= [1, 2, 4, 8, 13, 20, 30]
SPEEDS				= [0.73, 1.46, 2.19]
CONDITIONS			= list()


def FillConditionList():
	for angle in ANGLES:
		for frequency in FREQUENCIES:
			for speed in SPEEDS:
				CONDITIONS.append((angle, frequency, speed))

def PrintList(ml):
	for elt in ml:
		print(elt)

def LoadConditionTrials():
    cdt = pickle.load( open("conditionTrials.p", "rb") ) # read binary
    return cdt

def main():
    FillConditionList()
    #PrintList(CONDITIONS)
    conditionTrials = LoadConditionTrials()
    """
    print(len(conditionTrials))
    print("\n")
    for i in range(10):
        print(len(conditionTrials[i]))
        print(len(conditionTrials[i][0]))
        print(len(conditionTrials[i][0][0]))
        print("\n")

    """
    ci = 1
    for condition in conditionTrials:
        print("Condition " + str(ci))
        ct = 1
        for trial in condition:
            print("Trial " + str(ct))
            PrintList(trial)
            ct += 1
        ci += 1
    """
	points = np.random.rand(30, 2)
	hull = ConvexHull(points)
	#print(points)
	plt.plot(points[:,0], points[:,1], 'o')

	for simplex in hull.simplices:
		plt.plot(points[simplex, 0], points[simplex, 1], 'r-')


	print(hull.area)

	plt.show()
	"""

if __name__ == "__main__":
	main()
