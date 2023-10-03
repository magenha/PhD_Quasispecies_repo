"""
Iterate over clade abundance file and save abundance for each haplotype in the GN.
So that, we must extract abundance of each of the haplotypes by clade, evaluate
if that haplotype was already identified and (i) yes: add the clade abundance to its
GN abundance (ii) no: create new haplo entry and set its abundance.
Save a JSON with the GN_abundance dictionary: keys are nodeIDs and values abundances
"""

# Imports:
import numpy as np;
import os, sys;
import json

# Function to save JSON file given a dictionary
def saveJSON(path,fileName, dictionary):
    # Check whether the specified path exists or not. If not, create directory
    if not os.path.exists(path):
      os.makedirs(path)

    path_to_file = path + "/" +fileName
    outFile = open(path_to_file, 'w')
    json.dump(dictionary,outFile,indent=2,sort_keys=True)
    outFile.close()

# Load variable values
#proteinName = str(sys.argv[1])
#startAA = str(sys.argv[2])
#endAA = str(sys.argv[3])

# Defining InOut data path
dataPathInOut = "./Data/Abundances/";

# Loading fileNames in list
allFileNames = os.listdir(dataPathInOut)
# Init dictionary with GN abundances
GN_abundance = {}

# Loop over data files:
for thisFileName in allFileNames:
	print("Processing file: " + thisFileName);
	# Load next data file:
	fIn = open(os.path.join(dataPathInOut, thisFileName), 'r');
	dL = fIn.read().splitlines();
	fIn.close()
	for line in dL:
		nID = str(line.split(', ')[0])
		abundance = int(line.split(', ')[1])
		try:
			GN_abundance[nID] += abundance
		except KeyError:
			GN_abundance[nID] = abundance

# save dictionary with GN abundances
saveJSON(dataPathInOut,'haplotypes_GN.json',GN_abundance)
