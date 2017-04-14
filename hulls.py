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

def SliceTrials(condition, window):
	slices = list()
	for trial in condition:
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

def GetAllSlices(conditionTrials, window):
	allSlices = list()
	for condition in conditionTrials:
		allSlices.append(SliceTrials(condition, window))
	return(allSlices)

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

def PrintSlices(slices):
	cs = 0
	for mySlice in slices:
		print("Slice " + str(cs))
		PrintList(mySlice)
		cs += 1

def PrintAllSlices(allSlices):
	cd = 1
	for conditionSlices in allSlices:
		print("Condition " + str(cd))
		PrintSlices(conditionSlices)
		cd += 1

def PlotPoints(points):
	plt.plot(points[:,0], points[:,1], 'o')
	for simplex in hull.simplices:
		plt.plot(points[simplex, 0], points[simplex, 1], 'r-')
	plt.show()

def GetConditionHullAreas(condition):
	conditionAreas = list()
	sl = 1
	for mySlice in condition:
		points = mySlice[:,-2:]
		#print(points)
		hull = ConvexHull(points)
		conditionAreas.append(hull.area)
		sl += 1
	return(conditionAreas)



def CullSlices(allSlices):
	newSlices = list()
	for condition in allSlices:
		conditionSlices = list()
		for mySlice in condition:
			newX = False
			newY = False
			firstX = mySlice[0, -2]
			firstY = mySlice[0, -1]
			for line in mySlice:
				if firstX != line[-2]:
					newX = True
				if firstY != line[-1]:
					newY = True
			if newX and newY:
				conditionSlices.append(mySlice)
		newSlices.append(conditionSlices)
	return(newSlices)

def Test(condition):
	mySlice = condition[60]
	points = mySlice[:,-2:]
	print(points)
	hull = ConvexHull(points)
	print(hull.area)

def GetAllAreas(allSlices):
	allAreas = list()
	cd = 1
	for condition in allSlices:
		#print("Condition " + str(cd))
		allAreas.append(GetConditionHullAreas(condition))
		cd += 1
	return(allAreas)

def PrintAreas(allAreas):
	cd = 1
	for condition in allAreas:
		print("Condition " + str(cd))
		for area in condition:
			print(area)
		cd += 1

def main():
	FillConditionList()
	#PrintList(CONDITIONS)
	conditionTrials = LoadConditionTrials()
	allSlices = GetAllSlices(conditionTrials, 0.2)
	#PrintAllSlices(allSlices)
	#print(len(allSlices))
	allSlices = CullSlices(allSlices)
	#print(len(allSlices))

	#GetConditionHullAreas(allSlices[26])

	allAreas = GetAllAreas(allSlices)
	PrintAreas(allAreas)


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
