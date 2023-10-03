"""

	script_testAccess.py: 

		This file is just a quick test to make sure that all properly formated files can be accessed. At the same time, it
		makes sure that all genotypes are of the same size. 

"""

# Imports: 
import numpy as np; 
import matplotlib.pyplot as plt; 
from mpl_toolkits import mplot3d; 
import networkx as nx; 
import os, sys; 
from copy import copy; 
from scipy import stats; 


# Defining paths: 
dataPathIn = "/home/brigan/Desktop/Research_CNB/DeepSequencing/Data/CollapsedSequences_WB/"; 
codePath = "/home/brigan/Desktop/Research_CNB/DeepSequencing/CodeGrandNetwork/"; 

# Loading data file names: 
fIn = open(os.path.join(codePath, "fileNames_r3.txt"), 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 

# Extracting length of first genotype, which is taken as reference length: 
thisFileName = allFileNames[0]; 
fIn = open(os.path.join(dataPathIn, thisFileName+".fasta"), 'r'); 
newLine = fIn.readline(); 
newLine = fIn.readline().rstrip(); 
refLength = len(newLine); 

# Loop over data files: 
for thisFileName in allFileNames: 
	print("Processing file: " + thisFileName); 
	# Load next data file: 
	fIn = open(os.path.join(dataPathIn, thisFileName+".fasta"), 'r'); 

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
			abundance = newLine.split('-')[1]; 
		# If reading a sequence, check that length is correct: 
		else: 
			fReadMetadata = True; 
			thisLength = len(newLine); 
			if (thisLength != refLength): 
				print("ACHTUNG!! at iRead=" + str(iRead) + ": " + str(thisLength) + " is different from " + str(refLength) + "!!! "); 
				print(newLine); 


	



