#!/usr/bin/env python3
import random
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.spatial import ConvexHull


END_OF_TIMES	= 20
ANGLES			= range(15, 181, 15)
FREQUENCIES		= [1,2] + list(range(4, 64, 4))
#ANGLES			= [60]
#FREQUENCIES		= [8, 16]
SPEEDS			= [2.19]
CONDITIONS		= [(2.19, 0, 1)]

def FillConditionList():
	for speed in SPEEDS:
		for angle in ANGLES:
			for frequency in FREQUENCIES:
				CONDITIONS.append((speed, angle, frequency))

def PrintList(ml):
	for elt in ml:
		print(elt)

def RotateDirections2D(tDir, angleBound):
	dx, dy = tDir
	angleBound = math.radians(angleBound)
	angle = random.uniform(-angleBound, angleBound)
	#angle = angleBound
	#print("angle: ", math.degrees(angle))
	cosa = math.cos(angle)
	sina = math.sin(angle)
	tDir[0] = dx*cosa - dy*sina
	tDir[1] = dx*sina + dy*cosa
	return(tDir)


def MoveTargets2D(condition):
	speed, angle, frequency = condition
	pos = np.array((0.0, 0.0))
	tDir = np.array((1.0, 0.0))
	time = 0.0
	period = 1.0/frequency # not safe if nil frequency
	#if angle == 0.0:
		#return( np.array( [[0.0, 0.0, 0.0], [END_OF_TIMES, END_OF_TIMES*speed, 0]] ) )
	#print("period: ", period)
	deltaTime = 0.001
	real_eot = END_OF_TIMES
	nbLines = int(real_eot/deltaTime) #+ 1

	#print("nbLines: ", nbLines)
	positions = np.zeros((nbLines, 3))
	line = 1
	lastRotation = 0.0
	while time < real_eot and line < nbLines:
		time += deltaTime
		pos += tDir * deltaTime * speed
		positions[line] = np.append(time, pos)
		if time - lastRotation >= period:
			tDir = RotateDirections2D(tDir, angle)
			lastRotation = time
		line += 1
	return(positions)

def FileNameFromCondition(condition):
	speed, angle, frequency = condition
	fname = "trajsForManuscript/synTraj_219_" + str(angle) + "_" + str(frequency) + ".pdf"
	return(fname)

def PlotPoints(positions, condition):
	positions = np.array(positions)
	plt.axis("equal")		# Same scale on both axes, or screwed up perception of angles
	plt.plot(positions[:,0], positions[:,1], 'o')
	#for simplex in hull.simplices:
		#plt.plot(positions[simplex, 0], positions[simplex, 1], 'r-')
	#plt.show()
	plt.savefig(FileNameFromCondition(condition), bbox_inches="tight")
	plt.clf() # clears the plot so that I can create a new one from a clean basis

def main():
	FillConditionList()
	trajectories = list()
	cd = 1
	for condition in CONDITIONS:
		print("Processing condition ", cd, " out of ", len(CONDITIONS))
		traj = MoveTargets2D(condition)
		trajectories.append(traj)
		PlotPoints(traj[:,1:], condition)
		cd += 1
	pickle.dump(trajectories, open("trajectories_for_manuscript.p", "wb")) # write binary
	#PrintList(trajectories[0])
	#PrintList(positions)
	#hull = ConvexHull(positions[:,1:])



if __name__ == "__main__":
	main()
