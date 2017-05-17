#!/usr/bin/env python3
import random
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.spatial import ConvexHull


END_OF_TIMES	= 200
ANGLES			= [0, 30, 60, 90, 120, 180]
FREQUENCIES		= [1, 2, 4, 8, 13, 20, 30]
SPEEDS			= [0.73, 1.46, 2.19]
CONDITIONS		= list()
SMALLIFY		= False


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

"""
void PeriodicRotateDirections2D(float* directions) {
	float angle, dx, dy, cosa, sina;
	for(int i=0; i<CHOSEN_NB_TARGETS_2D; ++i) {
		// The angle should end up between range min and max.
		angle = RANGE_MIN + static_cast <float> (rand()) / (static_cast <float> ( RAND_MAX/(RANGE_MAX-RANGE_MIN) ) );
		dx = directions[3*i];
		dy = directions[3*i + 1];
		cosa = cos(angle);
		sina = sin(angle);
		directions[3*i]     = dx*cosa - dy*sina;
		directions[3*i + 1] = dx*sina + dy*cosa;
	}
}
"""

def MoveTargets2D(condition):
	#positions = list()

	speed, angle, frequency = condition
	pos = np.array((0.0, 0.0))
	tDir = np.array((1.0, 0.0))
	time = 0.0
	period = 1.0/frequency # not safe if nil frequency
	#if angle == 0.0:
		#return( np.array( [[0.0, 0.0, 0.0], [END_OF_TIMES, END_OF_TIMES*speed, 0]] ) )
	#print("period: ", period)
	deltaTime = 0.0001
	if SMALLIFY:
		deltaTime *= 10
		real_eot = END_OF_TIMES * 0.2
		nbLines = int(real_eot/deltaTime)
	else:
		real_eot = END_OF_TIMES
		nbLines = int(real_eot/deltaTime) #+ 1

	#print("nbLines: ", nbLines)
	positions = np.zeros((nbLines, 3))
	line = 1
	lastRotation = 0.0
	while time < real_eot and line < nbLines:
		time += deltaTime
		#print("tDir: ", tDir)
		#print("period: ", period)
		#print("speed: ", speed)
		pos += tDir * deltaTime * speed
		#print("pos: ", pos)
		positions[line] = np.append(time, pos)
		if time - lastRotation >= period:
			tDir = RotateDirections2D(tDir, angle)
			lastRotation = time
		#print("tDir and norm:", tDir, np.linalg.norm(tDir))
		line += 1
	return(positions)

def SetRanges(positions):
	minX = min(positions[:,1])
	maxX = max(positions[:,1])
	minY = min(positions[:,2])
	maxY = max(positions[:,2])
	rangeX = maxX - minX
	rangeY = maxY - minY
	xBigger = rangeX > rangeY
	diff = abs(rangeX-rangeY)
	if xBigger:
		minY -= diff/2
		maxY += diff/2
	else:
		minX -= diff/2
		maxX += diff/2
	return(minX, maxX, minY, maxY)

def PlotPoints(positions, hull):
	positions = np.array(positions)
	plt.axis("equal")		# Same scale on both axes, or screwed up perception of angles
	plt.plot(positions[:,0], positions[:,1], 'o')
	for simplex in hull.simplices:
		plt.plot(positions[simplex, 0], positions[simplex, 1], 'r-')
	plt.show()

def main():
	FillConditionList()
	trajectories = list()
	#positions = MoveTargets2D((0.6, 10, 13))
	#PrintList(positions)
	cd = 1
	for condition in CONDITIONS:
		print("Processing condition " +str(cd) + " out of " + str(len(CONDITIONS)))
		traj = MoveTargets2D(condition)
		trajectories.append(traj)
		cd += 1
	if SMALLIFY:
		pickle.dump(trajectories, open("mini_trajectories.p", "wb"))
	else:
		pickle.dump(trajectories, open("trajectories.p", "wb")) # write binary
	#PrintList(trajectories[0])
	#PrintList(positions)
	#hull = ConvexHull(positions[:,1:])
	#PlotPoints(positions[:,1:], hull)
	#print(hull.area)


if __name__ == "__main__":
	main()
