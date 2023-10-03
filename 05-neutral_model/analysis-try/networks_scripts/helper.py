"""

	helper.py: 

		Helper file to work with sequences. 

"""


import numpy as np; 
import matplotlib.pyplot as plt
import networkx as nx; 
import os, sys; 
from copy import copy; 
from scipy import stats; 


def findRelevantLoci(sequencesList): 
	"""	findRelevantLoci function: 

			This function runs over all sequences in sequencesList (all with same length) and finds relevant loci (i.e. loci in
			which there is at least one mutant). It returns the list of relevant loci as well as those where all sequences agree.
			This should make our lives easier when working with sequences. 

			Inputs: 
				>> sequencesList: From which we wish to find relevant loci. 

			Returns: 
				<< relevantLoci: Loci in which there is at leas one mutant. 
				<< irrelevantLoci: Loci in which there is no mutant in the whole sequence. 
					< ACHTUNG!! Not implemented yet. It's easier to generate if needed than to pass it along. 

	"""

	relevantLoci = []; 

	seqsLength = len(sequencesList[0]); 
	allLoci = [ii for ii in range(seqsLength)]; 
	unknownLoci = copy(allLoci); 

	seq1 = sequencesList[0]; 
	for seq2 in sequencesList[1::]: 
		iDiffer = np.squeeze(np.argwhere([n1!=n2 for (n1, n2) in zip(seq1, seq2)]), 1); 
		if (np.any(iDiffer)): 
			relevantLoci += list(iDiffer); 

	return np.unique(relevantLoci); 


def hammingDistance(str1, str2): 
	"""	hammingDistance function: 

			Computes the Hamming distance between two strings. 

			Inputs: 
				>> str1, str2: Strings from which the Hamming distance will be computed. 

			Returns: 
				<< hDist: Hamming distance between strings str1 y str2. 
	
	"""

	hDist = sum(cc1 != cc2 for (cc1, cc2) in zip(str1, str2)); 

	return hDist; 


def findSeqAbundances(allSequencesList, uniqueSequencesList=[], normed=False): 
	"""	findSeqAbundances function: 

			This function returns the abundances of each sequences in uniqueSequencesList found in allSequencesList. The list
			sequencesList is expected to contain unique elements only. If not provided ([] by default), the unique sequences list
			is built from allSequencesList. 

			Inputs: 
				>> allSequencesList: Lost of sequences of which we wish to obtain the abundances. 
				>> uniqie: Lis of unique sequences. 

			Returns: 
				<< seqAbundances: List with the abundances of each 

	"""

	# Have unique elements been obtained before? 
	if (uniqueSequencesList==[]): 
		uniqueSequencesList = np.unique(allSequencesList); 

	# Getting abundances of unique elements: 
	seqAbundances = []; 
	for seq in uniqueSequencesList: 
		seqAbundances += [allSequencesList.count(seq)]; 

	# Do we want this quantity normalized? 
	if (normed): 
		seqAbundances = np.divide(np.array(seqAbundances), sum(seqAbundances)); 

	return seqAbundances; 


def buildNetwork(sequencesList, hammingDistanceThreshold=1): 
	"""	buildNetwork function: 

			This function returns a networkx graph object in which nodes are sequences and edges appear for sequences which are
			not further appart than hammingDistance mutations. 

			Inputs: 
				>> sequencesList: List with all sequences from which to build the network. 
				>> hammingDistanceThreshold=1: By default, network is built for point mutations. 

			Returns: 
				<< G: Graph object with sequences as nodes and links connecting point mutations. 

	"""

	G = nx.Graph(); 
	G.add_nodes_from(sequencesList); 

	edgesList = []; 
	for (iSeq0, seq0) in enumerate(G.nodes): 
		for seq1 in list(G.nodes)[iSeq0+1::]: 
			if (hammingDistance(seq0, seq1)<=hammingDistanceThreshold): 
				edgesList += [(seq0, seq1)]; 

	G.add_edges_from(edgesList); 

	return G; 


def binarizeSequence(sequence): 
	"""	binarizeSequence function: 

			This functino turns the given sequence into a binary string. Therefore, each nucleotide is substituted by a
			4-component boolean vector that has a 0 everywhere except in the position corresponding to the nucleotide letter,
			taking: A=[1, 0, 0, 0], C=[0, 1, 0, 0], G=[0, 0, 1, 0], and T=[0, 0, 0, 1]. This scheme uses more bits than the
			minimum required, but this avoids introducing spurious correlations between the nucleotide types. 

			Inputs: 
				>> sequence: To be binirized. 

			Returns: 
				<< bSeq: binary representation of the input sequence according to the previous coding scheme. 

	"""

	# Dummy dictionary to speed up the process: 
	bDict = {}; 
	bDict["A"] = [1, 0, 0, 0]; 
	bDict["C"] = [0, 1, 0, 0]; 
	bDict["G"] = [0, 0, 1, 0]; 
	bDict["T"] = [0, 0, 0, 1]; 

	bSeq = []; 
	for nucleotide in sequence: 
		bSeq += bDict[nucleotide]; 

	return bSeq; 

def binarizeSequences(sequences): 
	"""	binarizeSequences function: 

			This function calls sequentially the previous function, binarizeSequence(), to generate the binarized version of all
			sequences in a sample. 

			Inputs: 
				>> sequences: To be binirized. 

			Returns: 
				<< bSeqs: Numpy matrix with the binirized version of the original sequences. 

	"""

	bSeqs = []; 
	for sequence in sequences: 
		bSeqs += [binarizeSequence(sequence)]; 

	return np.array(bSeqs).T; 
