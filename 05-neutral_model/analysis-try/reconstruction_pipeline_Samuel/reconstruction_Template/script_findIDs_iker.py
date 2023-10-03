"""

	script_findIDs.py:

		This script builds the dictionary of unique genotypes present in the Grand Network (GN). This is just a file
		containing all the unique sequences across all files of a region. The ID is just the order in which a unique string
		appears. This ID will be used later to produce more manageable networks.

"""

# Imports:
import numpy as np;
import os, sys;

# Defining paths and other auxiliary string variables:
dataPathIn = "./Data/";
dataPathOut = "./Data/";
codePath = "./";

# Loading data file names:
fIn = open(os.path.join(dataPathIn, "fileNames.csv"), 'r');
allFileNames = fIn.read().splitlines();
fIn.close();

# Loop over data files:
genotypeList = [];
genotypeDict = {};
for thisFileName in allFileNames:
	# if thisFileName in ['haplotypes_19A_319_321','haplotypes_19B_319_321','haplotypes_20A_319_321']:
	print("Processing file: " + thisFileName);
	# Load next data file:
	fIn = open(os.path.join(dataPathIn, thisFileName+".fasta"), 'r');
	# newLine = fIn.readline().rstrip();
	# print(newLine)

	# Flags to manage reading flow:
	fRead = True;
	fReadMetadata = True;
	iRead = 0;

	# While valid reads are returned:
	while(fRead):
		iRead += 1;
		# Check if next read is valid:
		newLine = fIn.readline().rstrip();
		if (not newLine):
			fRead = False;
			break;

		# If reading data, split the abundance out:
		if (fReadMetadata):
			fReadMetadata = False;
		# If reading a sequence, check that length is correct:
		else:
			fReadMetadata = True;
			genotypeList += [newLine];
			try:
				genotypeDict[newLine] += 1
			except KeyError: # new haplotype
				genotypeDict[newLine] = 1

# print(genotypeList)
# print(genotypeDict)

# Finding out unique sequences and writing them to file:
uniqueGenotypes = np.unique(genotypeList);
fOut = open(os.path.join(dataPathOut, "nodesID_old.csv"), 'w');
for seq in uniqueGenotypes:
	fOut.write(seq + '\n');
fOut.close();

# fOut2 = open(os.path.join(dataPathOut, "nodesID_iker.csv"), 'w');
fOut2 = open(os.path.join(dataPathOut, "nodesID.csv"), 'w');
for seq in genotypeDict.keys():
	fOut2.write(seq + '\n');
fOut2.close();
