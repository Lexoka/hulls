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
MODE			= "speeds"

def FillConditionList():
	print(ANGLES, FREQUENCIES)
	for speed in SPEEDS:
		for angle in ANGLES:
			for frequency in FREQUENCIES:
				CONDITIONS.append((speed, angle, frequency))

def AreaFillConditionList():
	global ANGLES
	global FREQUENCIES
	global CONDITIONS
	ANGLES		= [1, 5, 15, 30, 45, 60, 90, 120, 180]
	FREQUENCIES	= [1, 2, 4, 8, 16, 60, 120, 240]
	CONDITIONS	= list()
	FillConditionList()

def SpeedFillConditionList():
	global ANGLES
	global FREQUENCIES
	global CONDITIONS
	global SPEEDS
	ANGLES		= [120]
	FREQUENCIES	= [2]
	SPEEDS		= [2.19/8, 2.19/4, 2.19/2, 2.19, 2.19*2, 2.19*4]
	CONDITIONS	= list()
	FillConditionList()

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
	deltaTime = 0.001
	real_eot = END_OF_TIMES
	nbLines = int(real_eot/deltaTime) #+ 1

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
	speed		= str(speed)
	angle		= str(angle)
	frequency	= str(frequency)
	speed = speed.replace('.', '_') # removing the dot, so the remainder isn't treated as a file extension
	if MODE == "comp":
		fname = "trajsForManuscript/" + speed + "_gen/synTraj_" + speed + "_" + angle + "_" + frequency + ".pdf"
	elif MODE == "areas":
		fname = "trajsForManuscript/areas/areaTraj_" + speed + "_" + angle + "_" + frequency + ".pdf"
	elif MODE == "speeds":
		fname = "trajsForManuscript/speeds/spTraj_" + speed + "_" + angle + "_" + frequency + ".pdf"
	return(fname)

def GetBounds(trajectories):
	lxBound = float("inf")
	lyBound = float("inf")
	hxBound = -float("inf")
	hyBound = -float("inf")

	for traj in trajectories:
		if min(traj[:,1]) < lxBound:
			lxBound = min(traj[:,1])
		if max(traj[:,1]) > hxBound:
			hxBound = max(traj[:,1])

		if min(traj[:,2]) < lyBound:
			lyBound = min(traj[:,2])
		if max(traj[:,2]) > hyBound:
			hyBound = max(traj[:,2])
	print(lxBound, hxBound, lyBound, hyBound)
	margin = 2
	return(lxBound - margin, hxBound + margin, lyBound - margin, hyBound + margin)


def PlotPoints(positions, condition, bounds):
	positions = np.array(positions)
	if MODE != "speeds":
		plt.axis("equal")		# Same scale on both axes, or screwed up perception of angles
	plt.plot(positions[:,0], positions[:,1], 'o')
	if MODE == "areas":
		hull = ConvexHull(positions)
		cnt = 1
		for simplex in hull.simplices:
			if cnt == 1:
				lbl = "%.6G" % hull.area + " cm²" # ^ significant digits
				plt.plot(positions[simplex, 0], positions[simplex, 1], 'r-', label=lbl, linewidth = 5)
			else :
				plt.plot(positions[simplex, 0], positions[simplex, 1], 'r-', linewidth = 5) # It puts a label on every line otherwise
			plt.legend()
			cnt += 1
	if MODE == "speeds":
		plt.axis( [	min(bounds[0], bounds[2]),
		 			max(bounds[1], bounds[3]),
					min(bounds[0], bounds[2]),
					max(bounds[1], bounds[3]) ] ) # Very dirty, should detect bounds instead
		#plt.show()
	plt.savefig(FileNameFromCondition(condition), bbox_inches="tight")
	plt.clf() # clears the plot so that I can create a new one from a clean basis

def main():
	if MODE == "comp":
		FillConditionList()
	elif MODE == "areas":
		AreaFillConditionList()
	elif MODE == "speeds":
		SpeedFillConditionList()

	trajectories = list()
	cd = 1
	for condition in CONDITIONS:
		print("Processing condition ", cd, " out of ", len(CONDITIONS))
		traj = MoveTargets2D(condition)
		trajectories.append(traj)
		cd += 1
	cd = 1
	bounds = GetBounds(trajectories)
	for cd in range(len(trajectories)):
		print("Drawing condition ", cd, " out of ", len(CONDITIONS) - 1)
		PlotPoints(trajectories[cd][:,1:], CONDITIONS[cd], bounds)

	if MODE == "comp":
		pickle.dump(trajectories, open("trajectories_for_manuscript.p", "wb")) # write binary


if __name__ == "__main__":
	main()
