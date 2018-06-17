#!/usr/bin/env python3
import sys
import random
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.spatial import ConvexHull


END_OF_TIMES		= 10								# Duration of trajectories in seconds
ANGLES				= [15]
FREQUENCIES			= [16]
SPEEDS				= [2.19]						#
AUTOCORRELATIONS	= [0.0, 0.4, 0.6, 0.8, 0.9, 0.95, 0.975, 0.9875, 0.99375, 0.996875, 0.9984375, 1.0]
CONDITIONS			= []
DELTA_TIME			= 0.001

MODE				= "autocorr"				# Sets the mode. Different modes will generate different sets of pictures.
												# It's not very pretty and there probably should be a .config file and/or launch options instead.

INIT_ANGLE			= 0
_previous_angle		= INIT_ANGLE


# Simply builds a condition list based on supplied lists of S, A and F.
def FillConditionList():
	#print(ANGLES, FREQUENCIES)
	for speed in SPEEDS:
		for angle in ANGLES:
			for frequency in FREQUENCIES:
				for autocorr in AUTOCORRELATIONS:
					CONDITIONS.append((speed, angle, frequency, autocorr))

def PrintList(ml):
	for elt in ml:
		print(elt)

# Rotates the direction vector tDir by an angle sampled between -angleBound and +angleBound
def RotateDirections2D(tDir, angleBound, autocorr):
	global _previous_angle
	dx, dy = tDir
	angleBound = math.radians(angleBound)
	angleRand = random.uniform(-angleBound, angleBound)
	angle = _previous_angle*autocorr + angleRand*(1.0 - autocorr)
	_previous_angle = angle
	cosa = math.cos(angle)
	sina = math.sin(angle)
	tDir[0] = dx*cosa - dy*sina
	tDir[1] = dx*sina + dy*cosa
	return(tDir)

# Simply generates a trajectory, based on the parameters included in condition.
# The trajectory lasts END_OF_TIMES seconds, with a time step dictated by DELTA_TIME.
# Returns a numpy array with 3 columns: time, x, y; and a row for each time step.
def MoveTargets2D(condition):
	speed, angle, frequency, autocorr = condition
	pos = np.array((0.0, 0.0))
	tDir = np.array((1.0, 0.0))
	time = 0.0
	period = 1.0/frequency # not safe if nil frequency
	real_eot = END_OF_TIMES # not really useful here but whatever
	nbLines = int(real_eot/DELTA_TIME) #+ 1

	#time	x	y
	positions = np.zeros((nbLines, 3))
	line = 1
	lastRotation = 0.0
	while time < real_eot and line < nbLines:
		time += DELTA_TIME
		pos += tDir * DELTA_TIME * speed
		positions[line] = np.append(time, pos)
		if time - lastRotation >= period:
			tDir = RotateDirections2D(tDir, angle, autocorr)
			lastRotation = time
		line += 1
	return(positions)

# Pretty much a disposable function that generates convenient file names based on the condition and the mode.
# Maybe something tidier and more flexible would be useful, but I doubt it would be worth the effort.
def FileNameFromCondition(condition, filter=0, window=0, factor=0):
	speed, angle, frequency, autocorr = condition
	speed		= str(speed)
	angle		= str(angle)
	frequency	= str(frequency)
	autocorr	= str(autocorr)

	speed		= speed.replace('.', '_') # removing the dot, so the remainder isn't treated as a file extension
	autocorr	= autocorr.replace('.', '_')
	if MODE == "autocorr":
		fname = "autocorr/" + speed + "_" + MODE + "_" + speed + "_" + angle + "_" + frequency+ "_" + autocorr
	fname += ".pdf"
	return(fname)

# This returns the minimum and maximum Xs and Ys for all trajectories combined.
# It's useful if you want to plot all trajectories using the same scale.
def GetBounds(trajectories):
	lxBound = float("inf")	# Simple, safe init
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
	#print(lxBound, hxBound, lyBound, hyBound)
	margin = 0.5
	return(lxBound - margin, hxBound + margin, lyBound - margin, hyBound + margin)

# This generates the figures as nice little PDFs for a given trajectory (positions), a condition and bounds.
# Depending on the mode, the bounds may not be used. It's kind of dirty to require them anyway but it will do for now.
def PlotPoints(positions, condition, bounds, filter=0, window=0, factor=0):
	positions = np.array(positions)
	plt.axis("equal")		# Same scale on both axes, or screwed up perception of angles
	plt.plot(positions[:,0], positions[:,1], '-')
	#plt.axis( [	min(bounds[0], bounds[2]),		# I know this looks weird but it basically just enforces the same bounds on both axes.
	#			max(bounds[1], bounds[3]),		# If you don't do this, the trajectory can look heavily distorted.
	#			min(bounds[0], bounds[2]),		# Technically, if you wanted to completely avoid distortion, you'd have to generate square pictures.
	#			max(bounds[1], bounds[3]) ] )	# It's probably not hard to do with PyPlot but it looks OK as is.

	plt.savefig(FileNameFromCondition(condition, filter, window, factor), bbox_inches="tight")
	plt.clf() # clears the plot so that I can create a new one from a clean basis


def main():
	# First, set up the conditions as required by each mode, or rather by each purpose for which a mode is chosen.
	global _previous_angle
	if MODE == "autocorr":
		FillConditionList()

	# Empty list of trajectories for init.
	trajectories = []

	cd = 0
	# First, we build all the trajectories.
	for condition in CONDITIONS:
		print("Processing condition ", cd+1, " out of ", len(CONDITIONS))
		_previous_angle = INIT_ANGLE
		traj = MoveTargets2D(condition)
		trajectories.append(traj)
		cd += 1

	# Now, we get the trajectories' bounds. They may not be needed but it's simpler to get them anyway, plus it's very cheap.
	bounds = GetBounds(trajectories)

	# And now we make the pretty pictures and save them.
	for cd in range(len(trajectories)):
		print("Drawing condition ", cd+1, " out of ", len(CONDITIONS))
		PlotPoints(trajectories[cd][:,1:], CONDITIONS[cd], bounds)

	# OK so in some modes it can be useful to save the trajectories with pickle, because in extreme cases they may take a while to generate.
	# For the purposes of generating pictures it has yet to prove useful, but who knows?
	if MODE == "comp":
		pickle.dump(trajectories, open("trajectories_for_manuscript.p", "wb")) # write binary

# So it launches the main function if you call the script, but not if you include it as a module; which you probably wouldn't want to do, but you could.
if __name__ == "__main__":
	main()
