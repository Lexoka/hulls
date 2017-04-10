#!/usr/bin/env python3

import numpy as np
import pickle

ANGLES				= [0, 30, 60, 90, 120, 180]
FREQUENCIES			= [1, 2, 4, 8, 13, 20, 30]
SPEEDS				= [0.73, 1.46, 2.19]
CONDITION_TRIALS	= list()
CONDITIONS			= list()
FOLDERS				= list()
FILES				= list()

def FillConditionList():
	for angle in ANGLES:
		for frequency in FREQUENCIES:
			for speed in SPEEDS:
				CONDITIONS.append((angle, frequency, speed))
				CONDITION_TRIALS.append([])

def CondToString(cond):
	a, f, s = cond
	return(str(a) + "_" + str(f) + "_" + str(s))

def FillFileList():
	for folder in FOLDERS:
		#base = "all_logs/target_and_cursor/" + folder + "/subject_" + folder + "_"
		base = folder + "/subject_" + folder + "_"
		for condition in CONDITIONS:
			fname = base + CondToString(condition) + ".csv"
			FILES.append(fname)

def GoThroughFolders():
	for i in range(5,8):
		FOLDERS.append(str(i).zfill(2))

def PrintList(ml):
	for elt in ml:
		print(elt)

def ReadFile(filename):
	#print(filename)
	with open(filename, "r") as infile:
		content = infile.readlines()
		# First, we remove the first line, which only includes headers.
		content = content[1:]
		for i in range(len(content)):
			content[i] = content[i].split()
			#content[i] = content[i][:9]					# removing text messages at the end
			content[i] = content[i][0:2] + content[i][6:8]	# only keep trial + time + tX + tY
			content[i] = list(map(float, content[i]))
	npConverted = np.array(content)
	return(npConverted)

def AddTrials(fileNb, condlog):
	lastBound = 0
	currentEntry = 0
	for currentTrial in range(1,5):
		while currentEntry < len(condlog) and condlog[currentEntry][0] == currentTrial:
			currentEntry += 1
		CONDITION_TRIALS[fileNb].append(condlog[lastBound:currentEntry])
		lastBound = currentEntry

def PrintTrials(trials):
	nb = 1
	for trial in trials:
		print("Trial " + str(nb))
		print(trial[0])
		print(trial[-1])
		nb += 1

def BuildTrials():
	fileNb = 0	# also condition number
	for fname in FILES:
		# Getting the contents of the file as an np array
		conditionLog = ReadFile(fname)

		# Separating the 4 trials of each condition into separate np arrays,
		# and adding them to the right list in CONDITION_TRIALS
		AddTrials(fileNb, conditionLog)
		fileNb += 1

		# reset because we have one file per condition, but also per subject,
		# so when we reach len(CONDITIONS), that means we're done with the
		# current subject, on to the next one, and back to condition 0
		if fileNb == len(CONDITIONS):
			fileNb = 0

def WriteConditionTrials():
	with open("CONDITION_TRIALS.txt", "w") as outfile:
		for cond in CONDITION_TRIALS:
			for trial in cond:
				for line in trial:
					outfile.write(str(line) + "\n")

def SaveConditionTrials():
	pickle.dump(CONDITION_TRIALS, open("conditionTrials.p", "wb")) # write binary

def main():
	GoThroughFolders()
	FillConditionList()
	FillFileList()
	BuildTrials()
	SaveConditionTrials()

	"""
	for cdt in CONDITION_TRIALS:
		for tr in cdt:
			print(len(tr))
	"""

	#PrintList(CONDITIONS)
	#PrintList(FILES)
	#PrintList(ReadFile(FILES[0]))
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
