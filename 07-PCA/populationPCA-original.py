"""

	populationPCA.py: 

		This script loads the abundance that each condition presents in each of the nodes of a subnetwork. This array (one per
		experimental condition) is used to implement PCA and visualize how the population moves over genotype space. 

"""

# Imports: 
import numpy as np; 
import matplotlib.pyplot as plt; 
from mpl_toolkits import mplot3d; 
import networkx as nx; 
import os, sys; 

# For 3D scatter: 
from mpl_toolkits.mplot3d import Axes3D; 




## Preliminary commands: 

# Variables to assit on loading (e.g. number of nodes, etc): 
nLoad = 10; 
regionName = "r1"; 
dicFolder = {}; 
dicFolder["r1"] = "Reg1/"; 
dicFolder["r2"] = "Reg2/"; 
dicFolder["r3"] = "Reg3/"; 
fNormalize = True; 


# Defining paths and names: 
dataPathBase = "/home/brigan/Desktop/Research_CNB/DeepSequencing/Data/GrandNetwork/" + dicFolder[regionName]; 
abundacesPath = os.path.join(dataPathBase, "Abundances/"); 
figPathOut = "/home/brigan/Desktop/Research_CNB/DeepSequencing/CodeVisualizeGrandNetwork/Pics/PlotsPCA_"+str(regionName)+"_nNodes-"+str(nLoad)+'/'; 
strCall = "mkdir "+figPathOut; 
os.system(strCall); 


## Finding up top represented genotypes across all processed experimental conditions: 

# A list used to store all genotype IDs and retrieve only the subnetwork: 
allGenotypeIDs = []; 

# We need to read the top nLoad genotypes from each experimental condition that we wish to process. 
# Looping over those files: 
fIn = open("listFProcess.csv", 'r'); 
# fIn = open("listFProcess_wExtraCase.csv", 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 
nConditions = len(allFileNames); 

for thisFileName in allFileNames: 
	# Reading genotype IDs from the abundances folders: 
	genotypeIDs_and_abundances = np.loadtxt(os.path.join(abundacesPath, thisFileName+".csv"), delimiter=', ', dtype=int); 
	genotypeIDs = genotypeIDs_and_abundances[0:nLoad,0]; 
	allGenotypeIDs += list(genotypeIDs); 

# List of unique genotypes: 
allGenotypeIDs = np.unique(allGenotypeIDs); 
nGenotypes = len(allGenotypeIDs); 
# Creating a dictionary that links each genotype to its position in this list (handy later!): 
dictIGenotype = {}; 
for (iID, gID) in enumerate(allGenotypeIDs): 
	dictIGenotype[gID] = iID; 


## Loading the abundances of all top genotypes (including those which, within a condition, are not top): 

# Now we need to sort the abundances that each file presents in each of the nodes as a nConditions x nGenotypes matrix.
# Each column of this matrix will be an array of abundances to proceed with PCA. nGenotypes > nLoad, meaning that some
# of the entries might be zero. For simplicity, we initialize the corresponding matrix to be all zeros: 
abundanceArrays = np.zeros((nGenotypes, nConditions)); 

# Loading abundances for each experimental condition: 
for (iCondition, thisFileName) in enumerate(allFileNames): 
	fInName = os.path.join(abundacesPath, thisFileName+".csv"); 
	# Erasing the dummy file just in case: 
	strCall = "rm dummy.csv"; 
	os.system(strCall); 

	for (iGenotype, thisGenotype) in enumerate(allGenotypeIDs): 
		# Selecting the line in the abundances file that starts with this genotype ID: 
		strToGrep = "'^" + str(thisGenotype) + ", '"; # ACHTUNG!! Necessary to add the ", " to prevent grepping thisGenotype as a substring. 
		strCall = "grep -m 1 " + strToGrep + ' ' + fInName + " >> dummy.csv"; 
		os.system(strCall); 

	# Those genotypes represented with abundance different from zero are stored in the matrix: 
	genotypeIDs_and_abundances = np.loadtxt("dummy.csv", delimiter=', ', dtype=int); 
	for (iID, gID) in enumerate(genotypeIDs_and_abundances[:,0]): 
		abundanceArrays[dictIGenotype[gID], iCondition] = genotypeIDs_and_abundances[iID,1]; 


# Plotting the abundances in logarithmic scale: 
fig = plt.figure(); 
plt.imshow(np.log(abundanceArrays)); 
plt.colorbar(); 
fig.savefig(figPathOut + "colomapAbundances.pdf"); 



## Proceeding with the PCA analysis: 

# Do we want to look at normalized populations? 
if (fNormalize): 
	abundaceSum = np.sum(abundanceArrays, 0); 
	abundanceArrays = np.divide(abundanceArrays, np.repeat(np.array([abundaceSum]), nGenotypes, 0)); 

# Standardizing distro: 
abundanceMean = np.mean(abundanceArrays, 1); 
abundanceStd = np.std(abundanceArrays, 1); 
abundanceArrays = abundanceArrays - np.transpose(np.repeat(np.array([abundanceMean]), nConditions, 0)); 
abundanceArrays = np.divide(abundanceArrays, np.transpose(np.repeat(np.array([abundanceStd]), nConditions, 0))); 

# Computing correlation matrix and diagonalizing: 
abundancesCov = np.cov(abundanceArrays); 
(eigVals, eigVects) = np.linalg.eig(abundancesCov); 
eigVects = eigVects.real; 

'''
	MISSING: 
	-Needs to select the n-components biggest values
	-Needs to select the components associated with those values
'''

# Plotting eigenvalues matrix: 
fig = plt.figure(); 
plt.plot(eigVals); 
fig.savefig(figPathOut + "eigenvalues.pdf"); 

#Plotting covariance matrix: 
fig = plt.figure(); 
plt.imshow(eigVects); 
fig.savefig(figPathOut + "eigenvectors.pdf"); 

# Projecting data into abundance eigenspace:
abundanceArrays_ = np.dot(np.transpose(abundanceArrays), abundanceArrays); 
'''
	CONCEPT ERROR:
	-Multipllication by the components matrix only, to reduce the matrix.
'''

## # Quick visualization (improved below): 
# fig = plt.figure(); 
# plt.scatter(abundanceArrays_[0,:], abundanceArrays_[1,:]); 

# fig = plt.figure(); 
# ax = fig.add_subplot(111, projection='3d'); 
# ax.scatter(abundanceArrays_[0,:], abundanceArrays_[1,:], abundanceArrays_[2,:]); 


## Now we would like to make some cool projections of each of the conditions as they evolve over experimental time
## (i.e., over passages): 

# Grouping the indexes of each condition over their time course: 
# 	ACH! Remember that files go from more passages to less passages: !! 
conditionsList = ["30C", "33C", "37C", "40C", "43C"]; 
conditionsList_ = ["33C", "37C", "40C", "43C"]; 
dictIConditions = {}; 
dictIConditions["30C"] = [3, 2, 1, 0]; 
dictIConditions["33C"] = [7, 6, 5, 4]; 
dictIConditions["37C"] = [11, 10, 9, 8]; 
dictIConditions["40C"] = [15, 14, 13, 12]; 
dictIConditions["43C"] = [16, 17, 18]; # ACHTUNG!! Without extra case! 
# dictIConditions["43C"] = [17, 18, 19, 16]; # ACHTUNG!! With extra case! 

# Making a dictionary of colors for nice plot: 
dictStyleConditions = {}; 
dictStyleConditions["30C"] = 'kx-'; 
dictStyleConditions["33C"] = 'ro-'; 
dictStyleConditions["37C"] = 'bs-'; 
dictStyleConditions["40C"] = 'g+-'; 
dictStyleConditions["43C"] = 'y^-'; 

# Generating trajectories over PCA space: 
dictTrajectories1 = {}; 
dictTrajectories2 = {}; 
dictTrajectories3 = {}; 
dictTrajectories4 = {}; 
for thisCondition in conditionsList: 
	dictTrajectories1[thisCondition] = abundanceArrays_[0, dictIConditions[thisCondition]]; 
	dictTrajectories2[thisCondition] = abundanceArrays_[1, dictIConditions[thisCondition]]; 
	dictTrajectories3[thisCondition] = abundanceArrays_[2, dictIConditions[thisCondition]]; 
	dictTrajectories4[thisCondition] = abundanceArrays_[3, dictIConditions[thisCondition]]; 

# Generating 2D plots: 
fig = plt.figure(); 
for thisCondition in conditionsList: 
	plt.plot(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictStyleConditions[thisCondition]); 
plt.xlabel("PC1"); 
plt.ylabel("PC2"); 
fig.savefig(figPathOut + "trajectories2D.pdf"); 

fig = plt.figure(); 
for thisCondition in conditionsList_: 
	plt.plot(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictStyleConditions[thisCondition]); 
plt.xlabel("PC1"); 
plt.ylabel("PC2"); 
fig.savefig(figPathOut + "trajectories2D_wo30C.pdf"); 

# # Generating 2D plot: 
# fig = plt.figure(); 
# for thisCondition in conditionsList: 
# 	plt.plot(dictTrajectories2[thisCondition], dictTrajectories4[thisCondition], dictStyleConditions[thisCondition]); 
# plt.xlabel("PC2"); 
# plt.ylabel("PC4"); 
# fig.savefig(figPathOut + "trajectories2D_pc24.pdf"); 


# Generating 3D plot: 
fig = plt.figure(); 
ax = fig.add_subplot(111, projection='3d'); 
for thisCondition in conditionsList: 
	ax.plot3D(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictTrajectories3[thisCondition], dictStyleConditions[thisCondition]); 
ax.set_xlabel("PC1"); 
ax.set_ylabel("PC2"); 
ax.set_zlabel("PC3"); 
fig.savefig(figPathOut + "trajectories3D.pdf"); 

fig = plt.figure(); 
ax = fig.add_subplot(111, projection='3d'); 
for thisCondition in conditionsList_: 
	ax.plot3D(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictTrajectories3[thisCondition], dictStyleConditions[thisCondition]); 
ax.set_xlabel("PC1"); 
ax.set_ylabel("PC2"); 
ax.set_zlabel("PC3"); 
fig.savefig(figPathOut + "trajectories3D_wo30C.pdf"); 

# fig = plt.figure(); 
# ax = fig.add_subplot(111, projection='3d'); 
# for thisCondition in conditionsList: 
# 	ax.plot3D(dictTrajectories2[thisCondition], dictTrajectories3[thisCondition], dictTrajectories4[thisCondition], dictStyleConditions[thisCondition]); 
# ax.set_xlabel("PC2"); 
# ax.set_ylabel("PC3"); 
# ax.set_zlabel("PC4"); 
# fig.savefig(figPathOut + "trajectories3D_pc234.pdf"); 

plt.show(); 

