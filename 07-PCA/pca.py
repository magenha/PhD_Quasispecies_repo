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

add_log=''
if 0==1:
    add_log='_log'


# Defining paths and names: 
dataPathBase = "/home/samuel/Documents/PhD/Quasispecies/Data/"; 
#dataPathBase = "/home/samuel/Documents/PhD/Quasispecies/Data/";
abundacesPath = os.path.join(dataPathBase, "Sequences_filtered/N_rem_rem/"); 
#figPathOut = "./Pics-proteins/PlotsPCA_allFiles_"+str(regionName)+f"{add_log}_nNodes-"+str(nLoad)+'/'; 
figPathOut = "./Pics/PlotsPCA_allFiles_"+str(regionName)+f"{add_log}_nNodes-"+str(nLoad)+'/'; 
strCall = "mkdir "+figPathOut; 
os.system(strCall); 

def analyze_fname(a):
	'''
	Function that looks in file_name the     values of
	temperature and step
	
	Returns:  
		- temperature (Celsius degrees)
		- step
	'''
	if a[0] == 'c':
		step = int(a.split('-')[1][1:])
		t=43
	elif a[0] == '3':
		step = int(a.split('-')[2])
		t=30
	return t, step

## Finding up top represented genotypes across all processed experimental conditions: 

# A list used to store all genotype IDs and retrieve only the subnetwork: 
allGenotypeIDs = []; 

# We need to read the top nLoad genotypes from each experimental condition that we wish to process. 
# Looping over those files: 
fIn = open(os.path.join(dataPathBase + "fileNames.csv"), 'r'); 
allFileNames = fIn.read().splitlines(); 
fIn.close(); 
allFileNames.sort()
nConditions = len(allFileNames); 

index_30 = []
index_43 = []

for i,thisFileName in enumerate(allFileNames): 
	# Reading genotype IDs from the abundances folders: 
	genotypeIDs_and_abundances = np.loadtxt(os.path.join(abundacesPath, thisFileName+"._Nrem_rem.fasta"), delimiter=',', dtype=int); 
	genotypeIDs = genotypeIDs_and_abundances[0:nLoad,0]; 
	allGenotypeIDs += list(genotypeIDs); 
	t,step = analyze_fname(thisFileName)
	if t ==30:
		index_30.append(i)
	elif t==43:
		index_43.append(i)

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
	genotypeIDs_and_abundances = np.loadtxt("dummy.csv", delimiter=',', dtype=int); 
	for (iID, gID) in enumerate(genotypeIDs_and_abundances[:,0]): 
		if len(add_log):
			abundanceArrays[dictIGenotype[gID], iCondition] = np.log10(1.0+genotypeIDs_and_abundances[iID,1]); 
		else:
			abundanceArrays[dictIGenotype[gID], iCondition] = genotypeIDs_and_abundances[iID,1];


# Plotting the abundances in logarithmic scale: 
fig = plt.figure(); 
plt.imshow(np.log(abundanceArrays)); 
plt.colorbar(); 
fig.savefig(figPathOut + "colormapAbundances.pdf"); 

#Activate this to do PCA computes
if 1==1:
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

    # Plotting eigenvalues matrix: 
    #Scree plot
    fig = plt.figure(); 
    y = eigVals
    x = np.arange(1,len(list(eigVals))+1)
    plt.plot(x,y, 'o-'); 
    plt.ylabel('eigenvalue')
    plt.xlabel('component')
    plt.title('Scree plot')
    fig.savefig(figPathOut + "eigenvalues.pdf"); 

    #Plotting covariance matrix: 
    fig = plt.figure(); 
    plt.imshow(eigVects); 
    fig.savefig(figPathOut + "eigenvectors.pdf"); 

    # Projecting data into abundance eigenspace: 
    abundanceArrays_ = np.dot(np.transpose(abundanceArrays), abundanceArrays);
    #print('abundance arrays: ', abundanceArrays) 
    #print('abundance arrays_: ', abundanceArrays_)
    ## Now we would like to make some cool projections of each of the conditions as they evolve over experimental time
    ## (i.e., over passages): 


#Do plots
if 0==1:
    # Grouping the indexes of each condition over their time course: 
    # 	ACH! Remember that files go from more passages to less passages: !! 
    conditionsList = ["30C", "43C"]; 
    ancestorsList = ["ancestor2"]; 
    #conditionsList_ = ["30C", "33C", "37C", "40C", "43C", "43C_"]; 
    #ancestorsList = ["ancestor25", "ancestor2"]; 
    # conditionsList_ = ["33C", "37C", "40C", "43C"]; 
    
    dictIConditions = {}; 
    dictIConditions["30C"] = index_30;     #We need to find which columns of the matrix are these
    #dictIConditions["33C"] = [17, 16, 15, 14]; 
    #dictIConditions["37C"] = [21, 20, 19, 18]; 
    #dictIConditions["40C"] = [26, 25, 24, 23]; 
    dictIConditions["43C"] = index_43;         #We need to find which columns of the matrix are these
    #dictIConditions["G"] = [35, 36, 1]; 
    #dictIConditions["TS"] = [5, 4, 3, 2]; 
    #dictIConditions["43C_"] = [22, 30, 7]; 
    #dictIConditions["ancestor25"] = [11, 0, 6]; 
    dictIConditions["ancestor2"] = [index_30[0],index_43[1]];           #We need to find which columns of the matrix are these
    print('43',index_43)
    # Making a dictionary of colors for nice plot: 
    dictStyleConditions = {}; 
    dictStyleConditions["30C"] = 'bx-'; 
    #dictStyleConditions["33C"] = 'ro-'; 
    #dictStyleConditions["37C"] = 'bs-'; 
    #dictStyleConditions["40C"] = 'g+-'; 
    dictStyleConditions["43C"] = 'r^-'; 
    #dictStyleConditions["G"] = 'co-'; 
    #dictStyleConditions["TS"] = 'mo-'; 
    #dictStyleConditions["43C_"] = 'yo-'; 
    #dictStyleConditions["ancestor25"] = 'ko'; 
    dictStyleConditions["ancestor2"] = 'bo'; 

    # Generating trajectories over PCA space: 
    dictTrajectories1 = {}; 
    dictTrajectories2 = {}; 
    dictTrajectories3 = {}; 
    dictTrajectories4 = {}; 
    for thisCondition in conditionsList+ancestorsList: 
        dictTrajectories1[thisCondition] = abundanceArrays_[0, dictIConditions[thisCondition]]; 
        dictTrajectories2[thisCondition] = abundanceArrays_[1, dictIConditions[thisCondition]]; 
        dictTrajectories3[thisCondition] = abundanceArrays_[2, dictIConditions[thisCondition]]; 
        dictTrajectories4[thisCondition] = abundanceArrays_[3, dictIConditions[thisCondition]]; 

    # Generating 2D plots: 
    fig = plt.figure(); 
    plt.scatter(abundanceArrays_[0,:], abundanceArrays_[1,:]); 
    for thisCondition in conditionsList: 
        plt.plot(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictStyleConditions[thisCondition], label=f'{thisCondition}'); 
    for thisCondition in ancestorsList:
        plt.plot(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictStyleConditions[thisCondition], markersize=5, label=f'{thisCondition}'); 
    plt.plot(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictStyleConditions[thisCondition]); 
    plt.xlabel("PC1"); 
    plt.ylabel("PC2"); 
    plt.legend(loc='best')
    plt.tight_layout()
    fig.savefig(figPathOut + f"{nLoad}-trajectories2D.pdf"); 
    fig.savefig(figPathOut + f"{nLoad}-trajectories2D.png", format='png'); 

    

    fig = plt.figure(); 
    for thisCondition in conditionsList: 
        plt.plot(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictStyleConditions[thisCondition], label=f'{thisCondition}'); 
    for thisCondition in ancestorsList: 
        plt.plot(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictStyleConditions[thisCondition], markersize=5, label=f'{thisCondition}'); 
    plt.plot(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictStyleConditions[thisCondition]); 
    plt.xlabel("PC1"); 
    plt.ylabel("PC2"); 
    plt.legend(loc='best')
    plt.tight_layout()
    fig.savefig(figPathOut + f"{nLoad}-trajectories2D_classics.pdf");
    fig.savefig(figPathOut + f"{nLoad}-trajectories2D_classics.png", format='png'); 



    # Generating 3D plot: 
    fig = plt.figure(); 
    ax = fig.add_subplot(111, projection='3d'); 
    ax.scatter(abundanceArrays_[0,:], abundanceArrays_[1,:], abundanceArrays_[2,:]); 
    for thisCondition in conditionsList: 
        ax.plot3D(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictTrajectories3[thisCondition], dictStyleConditions[thisCondition], label=f'{thisCondition}'); 
    for thisCondition in ancestorsList: 
        ax.plot3D(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictTrajectories3[thisCondition], dictStyleConditions[thisCondition], markersize=5, label=f'{thisCondition}'); 
    ax.set_xlabel("PC1"); 
    ax.set_ylabel("PC2"); 
    ax.set_zlabel("PC3"); 
    ax.legend(loc='best')
    plt.tight_layout()
    fig.savefig(figPathOut + f"{nLoad}-trajectories3D.pdf"); 
    fig.savefig(figPathOut + f"{nLoad}-trajectories3D.png", format='png'); 

    fig = plt.figure(); 
    ax = fig.add_subplot(111, projection='3d'); 
    for thisCondition in conditionsList: 
        ax.plot3D(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictTrajectories3[thisCondition], dictStyleConditions[thisCondition], label=f'{thisCondition}'); 
    for thisCondition in ancestorsList: 
        ax.plot3D(dictTrajectories1[thisCondition], dictTrajectories2[thisCondition], dictTrajectories3[thisCondition], dictStyleConditions[thisCondition], markersize=5, label=f'{thisCondition}'); 
    ax.set_xlabel("PC1"); 
    ax.set_ylabel("PC2"); 
    ax.set_zlabel("PC3"); 
    ax.legend(loc='best')
    plt.tight_layout()
    fig.savefig(figPathOut + f"{nLoad}-trajectories3D_classics.pdf"); 
    fig.savefig(figPathOut + f"{nLoad}-trajectories3D_classics.png", format='png'); 

