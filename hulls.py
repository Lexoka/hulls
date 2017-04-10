#!/usr/bin/env python3

import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

ANGLES				= [0, 30, 60, 90, 120, 180]
FREQUENCIES			= [1, 2, 4, 8, 13, 20, 30]
SPEEDS				= [0.73, 1.46, 2.19]
CONDITIONS			= list()

"""
	conditionTrials:
		condition 0
			trial 0
				trial number				time		targetX	targetY
				1 to 4 (per subject basis)	x.x sec		0 to 1000, in pixels
			trial 1
			trial 2
			.
			.
			.
			trial 12 (4 per subject, 3 subjects at the moment)
		condition 1
		condition 2
		.
		.
		.
		condition N
"""

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

def SliceTrial(trial, window):
	slices = list()
	currentSlice = list()
	lastBeginning = trial[0,1]
	for line in trial:
		if line[1] - lastBeginning > window:
			currentSlice = np.array(currentSlice)
			slices.append(currentSlice)
			currentSlice = list()
			lastBeginning = line[1]
		currentSlice.append(line)
	return(slices)

def PrintConditionTrial(conditionTrials):
	ci = 1
	for condition in conditionTrials:
		print("Condition " + str(ci))
		ct = 1
		for trial in condition:
			print("Trial " + str(ct))
			PrintList(trial)
			ct += 1
		ci += 1

def main():
	FillConditionList()
	#PrintList(CONDITIONS)
	conditionTrials = LoadConditionTrials()
	"""
	res = SliceTrial(conditionTrials[0][0], 0.2)
	print(len(res))
	cs = 0
	for myslice in res:
		print("Slice " + str(cs))
		PrintList(myslice)
		cs += 1
	"""

	

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
