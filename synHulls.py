#!/usr/bin/env python3
import random
import math
import matplotlib.pyplot as plt
import numpy as np

END_OF_TIMES = 30

def PrintList(ml):
	for elt in ml:
		print(elt)

def RotateDirections2D(tDir, angleBound):
	dx, dy = tDir
	angleBound = math.radians(angleBound)
	angle = random.uniform(-angleBound, angleBound)
	angle = angleBound
	#print("angle: ", angle)
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
	positions = list()
	speed, angle, frequency = condition
	pos = np.array((0.0, 0.0))
	tDir = np.array((1.0, 0.0))
	time = 0.0
	if frequency == 0:
		period = END_OF_TIMES
	else:
		period = 1.0/frequency
	while time < END_OF_TIMES:
		time += period
		pos += tDir * period
		print(pos)
		positions.append(pos)
		tDir = RotateDirections2D(tDir, angle)
		#print(tDir, np.linalg.norm(tDir))
	return(positions)

def PlotPoints(positions):
	positions = np.array(positions)
	#plt.plot(positions[:,0], positions[:,1], 'o')
	plt.plot(positions)
	plt.show()

def main():
	positions = MoveTargets2D((0.81, 8.0, 2.0))
	PrintList(positions)
	PlotPoints(positions)


if __name__ == "__main__":
	main()
