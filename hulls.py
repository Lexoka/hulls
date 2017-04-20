#!/usr/bin/env python3

import numpy as np
import pickle
import matplotlib.pyplot as plt
import scipy.stats
from scipy.spatial import ConvexHull


ANGLES				= [0, 30, 60, 90, 120, 180]
FREQUENCIES			= [1, 2, 4, 8, 13, 20, 30]
SPEEDS				= [0.73, 1.46, 2.19]
CONDITIONS			= list()
WINDOWS				= [0.1, 0.15, 0.18, 0.19, 0.2, 0.21, 0.215, 0.22, 0.225, 0.23, 0.25, 0.3, 0.4, 0.6, 0.8]
WINDOWS				= [0.1, 0.2, 0.3, 0.4, 0.5]
WINDOWS				= [0.3, 0.7, 1.0]
DIAMETER			= 0.045245078
HDIAM				= DIAMETER/2.0

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
	for speed in SPEEDS:
		for angle in ANGLES:
			for frequency in FREQUENCIES:
				CONDITIONS.append((angle, frequency, speed))
				#print(angle, frequency, speed)

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
		print(mySlice[-1,1] - mySlice[0,1])
		PrintList(mySlice)
		cs += 1

def PrintAllSlices(allSlices):
	cd = 1
	for conditionSlices in allSlices:
		print("Condition " + str(cd))
		PrintSlices(conditionSlices)
		cd += 1

def PlotPoints(points, hull):
	plt.plot(points[:,0], points[:,1], 'o')
	for simplex in hull.simplices:
		plt.plot(points[simplex, 0], points[simplex, 1], 'r-')
	plt.show()

def SquarePoints(points):
	sqPoints = np.zeros((4*len(points), 2))
	i = 0
	for point in points:
		x,y = point
		lx = x - HDIAM
		rx = x + HDIAM
		by = y - HDIAM
		ty = y + HDIAM

		sqPoints[i]		= (lx, by)
		sqPoints[i+1]	= (lx, ty)
		sqPoints[i+2]	= (rx, by)
		sqPoints[i+3]	= (rx, ty)
		i += 4

	return(sqPoints)

def GetConditionHullAreas(condition, show, squarize = False):
	conditionAreas = list()
	sl = 1
	for mySlice in condition:
		points = mySlice[:,-2:]
		if squarize:
			points = SquarePoints(points)
		hull = ConvexHull(points)
		conditionAreas.append(hull.area)
		#print(hull.area)
		if show:
			PlotPoints(points, hull)

		sl += 1
	return(conditionAreas)

# Sometimes there are slices for which all values of X or Y may be equal, which
# boils down to a unidimensional space. This removes those slices.
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

# Get all areas, i.e., for all slices.
def GetAllAreas(allSlices):
	allAreas = list()
	cd = 1
	for condition in allSlices:
		#print("Condition " + str(cd))
		#allAreas.append(GetConditionHullAreas(condition, cd == 94))
		allAreas.append(GetConditionHullAreas(condition, False, False))
		cd += 1
	return(allAreas)

def GetAverageAreas(allAreas):
	averageAreas = list()
	averageStds = list()
	for condition in allAreas:
		average = np.mean(condition)
		std = np.std(condition)
		averageAreas.append(average)
		averageStds.append(std)
	return(averageAreas, averageStds)


def PrintAreas(allAreas):
	cd = 1
	for condition in allAreas:
		print("Condition " + str(cd))
		for area in condition:
			print(area)
		cd += 1

def WriteAverageAreas(averageAreas, window):
	fname = "areas_" + str(window) + ".txt"
	with open(fname, "w") as outfile:
		for area in averageAreas:
			outfile.write(str(area) + "\n")

def WriteAllAreas(allAreas, window):
	fname = "areas_all_" + str(window) + ".txt"
	with open(fname, "w") as outfile:
		cd = 1
		for condition in allAreas:
			outfile.write("Condition " + str(cd) + ":\n")
			cd += 1
			for area in condition:
				outfile.write(str(area) + "\n")


def ReadTimes():
	with open("average_results_all_speeds.csv", "r") as infile:
		contents = list()
		for line in infile:
			spLine = line.split(",")
			#spLine = line
			if spLine[0][0] != "#":
				contents.append(float(spLine[3]))
		#PrintList(contents)
	return(contents)

def main():

	selTimes = ReadTimes()
	#print(selTimes)
	FillConditionList()
	conditionTrials = LoadConditionTrials()

	lastThird = int(len(selTimes)/3)
	#print(lastThird)
	for window in WINDOWS:
		allSlices = GetAllSlices(conditionTrials, window)
		allSlices = CullSlices(allSlices)
		#PrintAllSlices(allSlices)
		allAreas = GetAllAreas(allSlices)
		averageAreas, averageStds = GetAverageAreas(allAreas)
		WriteAllAreas(allAreas, window)
		WriteAverageAreas(averageAreas, window)

		"""
		allSlices = CullSlices(allSlices)
		allAreas = GetAllAreas(allSlices)
		averageAreas, averageStds = GetAverageAreas(allAreas)
		WriteAverageAreas(averageAreas, window)
		print("Pearson for " + str(window))
		print(scipy.stats.pearsonr(selTimes, averageAreas))
		#print(scipy.stats.pearsonr(selTimes[-lastThird:], averageAreas[-lastThird:]))
		#print(scipy.stats.pearsonr(selTimes[lastThird:], averageAreas[lastThird:]))
		print("")
		"""



	"""
	points = np.random.rand(30, 2)
	print(points)
	print("")
	SquarePoints(points)
	"""

	"""
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
