#!/usr/bin/env python3
import sys
import random
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.spatial import ConvexHull


END_OF_TIMES	= 5								# Duration of trajectories in seconds
ANGLES			= [120]
FREQUENCIES		= [32]
SPEEDS			= [2.19]						#
CONDITIONS		= []
MODE			= "freqFilter"				# Sets the mode. Different modes will generate different sets of pictures.
											# It's not very pretty and there probably should be a .config file and/or launch options instead.

# Simply builds a condition list based on supplied lists of S, A and F.
def FillConditionList():
	#print(ANGLES, FREQUENCIES)
	for speed in SPEEDS:
		for angle in ANGLES:
			for frequency in FREQUENCIES:
				CONDITIONS.append((speed, angle, frequency))

def PrintList(ml):
	for elt in ml:
		print(elt)

# Rotates the direction vector tDir by an angle sampled between -angleBound and +angleBound
def RotateDirections2D(tDir, angleBound):
	dx, dy = tDir
	angleBound = math.radians(angleBound)
	angle = random.uniform(-angleBound, angleBound)
	cosa = math.cos(angle)
	sina = math.sin(angle)
	tDir[0] = dx*cosa - dy*sina
	tDir[1] = dx*sina + dy*cosa
	return(tDir)

# Simply generates a trajectory, based on the parameters included in condition.
# The trajectory lasts END_OF_TIMES seconds, with a time step dictated by deltaTime.
# Returns a numpy array with 3 columns: time, x, y; and a row for each time step.
def MoveTargets2D(condition):
	speed, angle, frequency = condition
	pos = np.array((0.0, 0.0))
	tDir = np.array((1.0, 0.0))
	time = 0.0
	period = 1.0/frequency # not safe if nil frequency
	deltaTime = 0.001
	real_eot = END_OF_TIMES # not really useful here but whatever
	nbLines = int(real_eot/deltaTime) #+ 1

	#time	x	y
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


def FilterTrajectory(traj, freq):
	if freq <= 0:
		return()
	period = 1.0/freq
	print(period)
	nbl, nbc = traj.shape
	print(nbl, nbc)
	#filtPositions = np.zeros((nbl,nbc))
	filtPositions = []
	lastStep = 0
	currTime = 0
	for i in range(nbl):
		currTime = traj[i,0]
		if (currTime - lastStep) >= period:
			#filtPositions[i] = traj[i]
			filtPositions.append(traj[i])
			lastStep = currTime
	filtPositions = np.array(filtPositions)
	print(len(filtPositions))
	return(filtPositions)

	#for line in filtPositions:
	#	print(line)

def GetDistance(traj):
	linDistancesCovered = np.diff(traj, axis=0)
	nbl, nbc = linDistancesCovered.shape
	twoDdist = np.zeros((nbl))
	#twoDdist = np.sqrt( linDistancesCovered[:,1]**2 + linDistancesCovered[:,2]**2 )
	twoDdist = np.sqrt( linDistancesCovered[:,0]**2 + linDistancesCovered[:,1]**2 ) # For PlotPoints, because it removes the time stamps
	#print(twoDdist)
	#print(twoDdist.sum())
	return(twoDdist.sum())

# Pretty much a disposable function that generates convenient file names based on the condition and the mode.
# Maybe something tidier and more flexible would be useful, but I doubt it would be worth the effort.
def FileNameFromCondition(condition, filter=0):
	speed, angle, frequency = condition
	speed		= str(speed)
	angle		= str(angle)
	frequency	= str(frequency)
	speed = speed.replace('.', '_') # removing the dot, so the remainder isn't treated as a file extension
	if MODE == "freqFilter":
		fname = "filteredTrajsChapFive/" + speed + "_" + MODE + "_" + speed + "_" + angle + "_" + frequency
		if filter > 0:
			fname += "_filter_" + str(filter)
		fname += ".pdf"
	elif MODE == "areas":
		fname = "trajsForManuscript/areas/areaTraj_" + speed + "_" + angle + "_" + frequency + ".pdf"
	elif MODE == "speeds":
		fname = "trajsForManuscript/speeds/spTraj_" + speed + "_" + angle + "_" + frequency + ".pdf"
	elif MODE == "autocorr":
		fname = "trajsForManuscript/autocorr/ac_" + speed + "_" + angle + "_" + frequency + ".pdf"
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
	margin = 1
	return(lxBound - margin, hxBound + margin, lyBound - margin, hyBound + margin)

# This generates the figures as nice little PDFs for a given trajectory (positions), a condition and bounds.
# Depending on the mode, the bounds may not be used. It's kind of dirty to require them anyway but it will do for now.
def PlotPoints(positions, condition, bounds, filter=0):
	positions = np.array(positions)
	plt.axis("equal")		# Same scale on both axes, or screwed up perception of angles
	distance = GetDistance(positions)
	lbl = "%.6G" % distance + " cm"
	plt.plot(positions[:,0], positions[:,1], '-', label=lbl)
	if MODE == "freqFilter":
		hull = ConvexHull(positions)
		cnt = 1
		for simplex in hull.simplices:
			if cnt == 1:
				lbl = "%.6G" % hull.area + " cm²" # 6 significant digits

				plt.plot(positions[simplex, 0], positions[simplex, 1], 'r-', label=lbl, linewidth = 5)
			else :

				plt.plot(positions[simplex, 0], positions[simplex, 1], 'r-', linewidth = 5) # It puts a label on every line otherwise,
																							# so I only set a label on the first one.
			plt.legend()
			cnt += 1
	if MODE == "speeds":
		plt.axis( [	min(bounds[0], bounds[2]),
		 			max(bounds[1], bounds[3]),
					min(bounds[0], bounds[2]),
					max(bounds[1], bounds[3]) ] )	# I know this looks weird but it basically just enforces the same bounds on both axes.
													# If you don't do this, the trajectory can look heavily distorted.
													# Technically, if you wanted to completely avoid distortion, you'd have to generate square pictures.
													# It's probably not hard to do with PyPlot but it looks OK as is.
		#plt.show()	# To display the picture in a sort of interactive window.
	plt.savefig(FileNameFromCondition(condition, filter), bbox_inches="tight")
	plt.clf() # clears the plot so that I can create a new one from a clean basis


def SaveApprox(allApprox, MY_ANGLES):
	with open("linApprox.csv", "w") as outfile:
		outfile.write("Angle	Frequency	m	p\n")
		for line in allApprox:
			a,f,m,p = line
			if a == MY_ANGLES[0]:
				outfile.write("\n") # Blank line for gnuplot with pm3d
			outfile.write(str(a) + "\t" + str(f) + "\t" + str(m) + "\t" + str(p) + "\n")


def main():
	FREQS_FOR_FILTER = [4,8,16]
	# First, set up the conditions as required by each mode, or rather by each purpose for which a mode is chosen.
	if MODE == "freqFilter":
		FillConditionList()
	else:
		pass

	# Empty list of trajectories for init.
	trajectories = []
	filtTrajectories = []
	cd = 0
	# First, we build all the trajectories.
	for condition in CONDITIONS:
		print("Processing condition ", cd, " out of ", len(CONDITIONS) - 1)
		traj = MoveTargets2D(condition)
		trajectories.append(traj)
		for FREQ_FOR_FILTER in FREQS_FOR_FILTER:
			filtTraj = FilterTrajectory(traj, FREQ_FOR_FILTER)
			filtTrajectories.append(filtTraj)
		cd += 1

	# Now, we get the trajectories' bounds. They may not be needed but it's simpler to get them anyway, plus it's very cheap.
	bounds = GetBounds(trajectories)

	# And now we make the pretty pictures and save them.
	for cd in range(len(trajectories)):
		print("Drawing condition ", cd, " out of ", len(CONDITIONS) - 1)
		PlotPoints(trajectories[cd][:,1:], CONDITIONS[cd], bounds)
		filtTrajIndex = 0
		for FREQ_FOR_FILTER in FREQS_FOR_FILTER:
			#PlotPoints(filtTrajectories[cd][:,1:], CONDITIONS[cd], bounds, FREQ_FOR_FILTER)
			PlotPoints(filtTrajectories[filtTrajIndex][:,1:], CONDITIONS[cd], bounds, FREQ_FOR_FILTER)
			filtTrajIndex += 1

	# OK so in some modes it can be useful to save the trajectories with pickle, because in extreme cases they may take a while to generate.
	# For the purposes of generating pictures it has yet to prove useful, but who knows?
	if MODE == "comp":
		pickle.dump(trajectories, open("trajectories_for_manuscript.p", "wb")) # write binary

# So it launches the main function if you call the script, but not if you include it as a module; which you probably wouldn't want to do, but you could.
if __name__ == "__main__":
	main()
