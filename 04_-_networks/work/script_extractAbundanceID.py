"""

	script_extractAbundanceID.py: 

		This script extracts the abundance and ID of each genotype in an experimental condition, and stores them in a file in
		which both ID and abundance (now in this order) become easy to access. 

		ACHTUNG!! Modifications underway to store normalized abundance as well. 

"""

# Imports: 
import numpy as np; 
import os, sys; 


# Defining paths and other auxiliary string variables: 
dataPathIn = os.getcwd() + '/Data/'
dataPathOut = os.getcwd() + '/Data/'
abundancesPathOut = dataPathOut + "Abundances/"
if (not os.path.exists(abundancesPathOut)):
	os.mkdir(abundancesPathOut)	


# Loading data file names: 
fIn = open(os.path.join(dataPathIn, "fileNames.csv"), 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 

# Setting name of ID file: 
fIDName = os.path.join(dataPathOut, "nodesID.csv"); 


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
	toWriteID = []; 
	toWriteAbundance = []; 
	Z1 = 0.; 
	Z2 = 0.; 
	while(fRead): 
		iRead += 1; 
		# Check if next read is valid: 
		newLine = fIn.readline().rstrip(); 
		if (not newLine): 
			fRead = False; 
			break; 

		# If reading data, split the abundance out: 
		if (fReadMetadata): 
			toWriteAbundance += [int(newLine.split('-')[1])]; 
			Z1 += toWriteAbundance[-1]; 
			Z2 += np.power(toWriteAbundance[-1], 2); 
			fReadMetadata = False; 

		# If reading a sequence, find out its ID in the grand network: 
		else: 
			fReadMetadata = True; 
			thisGenotype = newLine; 
			strOut = "grep -n " + thisGenotype + ' ' + fIDName + " > ./temp.csv"; 
			os.system(strOut); 
			fIn_ = open("./temp.csv", 'r'); 
			tempLine = fIn_.readline(); 
			fIn_.close(); 
			toWriteID += [int(tempLine.split(':')[0])-1]; 

	# Normalizing abundances: 
	# 	ACHTUNG!! It will be useful to compute two different normalizations: 
	# 		> In the first one, components add up to 1, thus this is a probability distribution. 
	# 		> In the second one, components are normalized to 1 with the Euclidean norm. 
	toWriteFraction = np.divide(toWriteAbundance, Z1); 
	toWriteComponent = np.divide(toWriteAbundance, np.sqrt(Z2)); 

	# File to save output: 
	fOut = open(os.path.join(dataPathOut, "Abundances", thisFileName + ".csv"), 'w'); 

	# Loop over data: 
	for (thisID, thisAbundance, thisFraction, thisComponent) in zip(toWriteID, toWriteAbundance, toWriteFraction, toWriteComponent): 
		strOut = str(thisID) + ', ' + str(thisAbundance) + ', ' + str(thisFraction) + ', ' + str(thisComponent) + '\n'; 
		fOut.write(strOut); 

	fOut.close(); 





	



