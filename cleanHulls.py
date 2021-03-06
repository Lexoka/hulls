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
#WINDOWS				= [0.5, 1.0]
WINDOWS				= [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
#WINDOWS				= [1.0]
DIAMETER			= 0.045245078
HDIAM				= DIAMETER/2.0
CUTOFF				= 4
CUTTING				= False
SMALLIFY			= False
MAKE_SQUARE_POINTS	= True


def FillConditionList():
	for speed in SPEEDS:
		for angle in ANGLES:
			for frequency in FREQUENCIES:
				CONDITIONS.append((speed, angle, frequency))

def PrintList(ml):
	for elt in ml:
		print(elt)

def SliceCondition(condition, window):
	slices = list()
	lastBeginning = condition[0,0]
	currentSlice = list()
	for line in condition:
		if line[0] - lastBeginning > window:
			currentSlice = np.array(currentSlice)
			slices.append(currentSlice)
			currentSlice = list()
			lastBeginning = line[0]
		currentSlice.append(line)
	return(slices)

def GetAllSlices(trajectories, window):
	allSlices = list()
	cd = 1
	for condition in trajectories:
		print("GetAllSlices, slicing condition " + str(cd) + "...")
		allSlices.append(SliceCondition(condition, window))
		print("GetAllSlices, done slicing condition " + str(cd) + ".")
		cd += 1
		if CUTTING and cd == CUTOFF:
			return(allSlices)
	return(allSlices)

def PlotPoints(positions, hull):
	positions = np.array(positions)
	#plt.axis("equal")		# Same scale on both axes, or screwed up perception of angles
	plt.plot(positions[:,0], positions[:,1], 'o')
	for simplex in hull.simplices:
		plt.plot(positions[simplex, 0], positions[simplex, 1], 'r-')
	plt.show()

def SquareTimedPoints(points):
	sqPoints = np.zeros((4*len(points), 3))
	i = 0
	for point in points:
		t,x,y = point
		lx = x - HDIAM
		rx = x + HDIAM
		by = y - HDIAM
		ty = y + HDIAM
		sqPoints[i]		= (t, lx, by)
		sqPoints[i+1]	= (t, lx, ty)
		sqPoints[i+2]	= (t, rx, by)
		sqPoints[i+3]	= (t, rx, ty)
		i += 4
	return(sqPoints)

def SquarePoints(points):
	#print(points)
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

def LinearTimedTrajectory(points):
	#print(points)
	if len(points) > 4:
		A = points[0]
		B = points[16]
		C = points[128]

		#print(A,B,C)

		At, Ax, Ay = A
		Bt, Bx, By = B
		Ct, Cx, Cy = C

		return(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By) < 0.0001)

def LinearTrajectory(points):
	#print(points)
	A = points[0]
	B = points[int(len(points)/2)]
	C = points[len(points) - 1]
	#print(A,B,C)

	Ax, Ay = A
	Bx, By = B
	Cx, Cy = C

	return(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By) < 0.0001)

def GetConditionHullAreas(condition, conditionParameters, show, squarize = False):
	speed, angle, frequency = conditionParameters
	conditionAreas = list()
	sl = 1
	for mySlice in condition:
		points = mySlice[:,-2:]
		if not squarize and angle == 0:
			conditionAreas.append(0.0)
		else:
			isLinear = (angle == 0 or LinearTrajectory(points))
			#print(angle == 0)
			#print(LinearTrajectory(points))
			#print(isLinear)
			if squarize:
				points = SquarePoints(points)
			if isLinear:
				#PrintList(points)
				bottomLeft = points[0]
				btlX, btlY = bottomLeft

				topRight = points[-1]
				tprX, tprY = topRight

				width = abs(tprX - btlX)
				height = abs(tprY - btlY)
				conditionAreas.append(width * height)
			else:
				#print(conditionParameters)
				#PrintList(points)
				hull = ConvexHull(points)
				conditionAreas.append(hull.area)
				if show:
					PlotPoints(points, hull)
		sl += 1
	return(conditionAreas)

def GetAllAreas(allSlices):
	allAreas = list()
	cd = 1
	for condition in allSlices:
		print("GetAllAreas, computing areas for condition " + str(cd) + "...")
		conditionParameters = CONDITIONS[cd-1]
		allAreas.append(GetConditionHullAreas(condition, conditionParameters, False, MAKE_SQUARE_POINTS))
		print("GetAllAreas, done computing areas for condition " + str(cd) + ".")
		cd += 1
		if CUTTING and cd == CUTOFF:
			return(allAreas)
	return(allAreas)

def WriteAverageAreas(averageAreas, window):
	fname = "clean_areas_" + str(window) + ".txt"
	with open(fname, "w") as outfile:
		for area in averageAreas:
			outfile.write(str(area) + "\n")

def WriteAllAreas(allAreas, window):
	fname = "clean_areas_all_" + str(window) + ".txt"
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
	return(contents)

def GetAverageAreas(allAreas):
	averageAreas = list()
	averageStds = list()
	cd = 1
	for condition in allAreas:
		print("GetAverageAreas, averaging for condition " + str(cd) + "...")
		average = np.mean(condition)
		std = np.std(condition)
		averageAreas.append(average)
		averageStds.append(std)
		cd += 1
	return(averageAreas, averageStds)

def main():
	selTimes = ReadTimes()
	FillConditionList()
	if SMALLIFY:
		fname = "mini_trajectories.p"
	else:
		fname = "trajectories.p"
	print("Loading", fname, "...")
	trajectories = pickle.load( open(fname, "rb") ) # read binary
	print(fname, "loaded.")

	"""
	for i in [0, 18]:
		traj = SquarePoints(trajectories[i])
		hull = ConvexHull(traj[:,1:])
		PlotPoints(traj[:,1:], hull)
		if LinearTrajectory(traj):
			print(42.42)
		else:
			print(hull.area)
	"""

	#"""
	for window in WINDOWS:
		allSlices = GetAllSlices(trajectories, window)
		allAreas = GetAllAreas(allSlices)
		print("Computing average areas...")
		averageAreas, averageStds = GetAverageAreas(allAreas)
		#WriteAllAreas(allAreas, window)
		print("Writing average areas...")
		WriteAverageAreas(averageAreas, window)
	#"""

if __name__ == "__main__":
	main()
