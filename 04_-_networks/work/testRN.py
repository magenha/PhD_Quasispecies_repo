"""

	testRN.py: 

		This is a script to quickly test that the renormalization has proceeded correctly. 

"""


def loadNodeCatalog(fInName): 
	"""	loadNodeCatalog function: 

			This function loads a node catalog. 

			Inputs: 
				>> fInName: File from where to load the node catalog. 

			Returns: 
				<< nodeCatalog: indicating which original sequences have been coarse grained into which nodes. 

	"""

	nodeCatalog = []; 
	allData = np.loadtxt(fInName).astype(int); 
	iRead = 0; 
	nNodes = allData[iRead]; 
	iRead += 1; 
	for ii in range(nNodes): 
		thisNSeqs = allData[iRead]; 
		iRead += 1; 
		theseSeqs = []; 
		for jj in range(thisNSeqs): 
			theseSeqs += [allData[iRead]]; 
			iRead += 1; 
		nodeCatalog += [theseSeqs]; 

	return nodeCatalog; 

def loadSequence(fInName, iSeq): 
	"""	loadSequence function: 

			This function loads the iSeq-th sequence from the file that contains all sequences. 

			Inputs: 
				>> fInName: From where to load sequences. 
				>> iSeq: Sequence number that we wish to load. 

			Returns: 
				<< seq: Sequence on the iSeq+1 position of file fInName. 

	"""

	strCall = "sed '" + str(iSeq+1) + "!d' "  + fInName + " > temp.txt"; 
	os.system(strCall); 

	fIn = open("temp.txt", 'r'); 
	seq = fIn.read().splitlines()[0]; 

	return seq; 



# Imports: 
import numpy as np; 
import os, sys; 

# Defining paths and other auxiliary string variables: 
allPath = "/home/brigan/Desktop/Research_CNB/Bio/DeepSequencing/Data/GrandNetwork/Toy2/"; 
pathGN = allPath + "GN/"; 
pathRN = allPath + "RN/"; 
fAllSeqs = allPath + "nodesID_toy2.csv"; 
fRepSeqs = "/home/brigan/Desktop/Research_CNB/Bio/DeepSequencing/Data/GrandNetwork/Toy2_RN/nodesID_toy2_RN.csv"; 

fCatalogName = pathRN + "nodeCatalog.csv"; 

nodeCatalog = loadNodeCatalog(fCatalogName); 
iIgnore = [(ii**2) for ii in range(10)]; 
iCheck = np.array([ii for ii in range(271) if ii not in iIgnore]); 
print(iIgnore); 
print(iCheck); 

# Now we have to check that all the nodes that have been collapsed into the same node present, indeed, the same sequence: 
fOut = open(fRepSeqs, 'w'); 
for node in nodeCatalog: 
	# Proceed only for nodes that coarse grain more than one sequence: 
	repSeq = loadSequence(fAllSeqs, node[0]); 
	for ii in iCheck: 
		fOut.write(repSeq[ii]); 
	fOut.write('\n'); 
	# if (len(node)>1): 
	# 	# Load representative sequence: 
	# 	for iSeq in np.arange(1, len(node)): 
	# 		thisSeq = loadSequence(fAllSeqs, node[iSeq]); 
	# 		if (not(repSeq[iCheck]==thisSeq[iCheck])): 
	# 			print("ACHTUNG!! "); 
	# 			print(node[0]); 
	# 			print('\t' + repSeq); 
	# 			print(node[iSeq]); 
	# 			print('\t' + thisSeq); 

fOut.close(); 


