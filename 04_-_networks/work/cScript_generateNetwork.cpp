/* 

	cScript_generateNetwork.cpp: 

		This is the central script to generate networks from the sequencing data. We have tried different ways to build these
		networks from the most exhausting (comparing all sequences nucleotide by nucleotide) to the nimbler (which uses a
		reference sequence and we store in which nucleotides each other differs, and those differences are used to guide
		comparison between other sequences). Each of this implementations is saved in individual files in the folder
		~/AlternativeImplementations/. 

		Now, in this very file we find the latest, working implementation -- which is usually the fastest one. 

		This script uses all pieces of information, ending up with the detailed divergences from a reference sequence, to
		reconstruct the genotype network. More specifically: 
			i) It uses the distance to a first reference sequence: 
				i.a) Said reference sequence will be an abundant sequence across all files (this has been found externally). 
					- For region r1, ID: 95217. 
					- For region r2, ID: 310589. 
					- For region r3, ID: 136622. 
					- ACHTUNG!! These IDs might change each time that the Grand Network reconstruction software is run. 
				i.b) For this reference, the distribution of sequences at a given distance is built. 
				i.c) We only need to compare sequences that differ at most in one more/less nucleotide than each other. 
			ii) It also uses the distance to a second reference sequence: 
				ii.a) This sequence is chosen as the one furthest away from the first reference. 
				ii.b) We look at the distances of sequences to this second reference to check whether they differ in more than one. 
			iii) It uses the specific locations at which each sequence differs from the first reference sequence. 
				- These locations are stored in a file for fast retrieval in case they are needed later. 

		ACHTUNG!! 
			>> This implementation will be improved to include annotated information about in which position does each point
			mutation separating two connected nodes take place. 
				> I have already done this and it turns out that we do not need it for renormalizing networks. 
			>> A pipeline will be build to run all the code from sequence, through finding abundances, to Grand Network. 

*/ 


// Imports: 
#include <iostream> 
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib> 
#include <vector>
#include <sstream>
#include "gNets.h"
using namespace std; 

#include <chrono> 
using namespace std::chrono; 


int main(int argc, char* argv[]){
	/* main routine: */ 

	// Defining variables and loading paths: 
	cout << "Defining variables and loading paths: " << endl; 
	// int refID=95217; 	// For r1! 
	// int refID=310589; 	// For r2! 
	// int refID=136622; 	// For r3! 
	int refID=0; 
	string netName="toy2"; 
	ostringstream dataPathIn, dataPathOut, sysStr, fInName, fOutName;
	dataPathIn.str(""); 
	dataPathIn << "/home/brigan/Desktop/Research_CNB/Bio/DeepSequencing/Data/GrandNetwork/Toy2/"; 
	gN::GNet gNet; 
	fInName.str(""); 
	fInName << dataPathIn.str() << "metadata_" << netName << ".csv"; 
	gNet.loadMetadata(fInName.str()); 
	cout << "Variables defined and paths loaded! " << endl; 

	// Creating out-paths: 
	cout << "Creating out-paths: " << endl; 
	dataPathOut.str(""); 
	dataPathOut << dataPathIn.str() << "GN/"; 
	sysStr.str(""); 
	sysStr << "rm -r " << dataPathOut.str(); 
	system(sysStr.str().c_str()); 
	sysStr.str(""); 
	sysStr << "mkdir " << dataPathOut.str(); 
	system(sysStr.str().c_str()); 
	cout << "Out-paths created! " << endl; 


	// Start chrono: 
	auto start = high_resolution_clock::now(); 

	// Finding auxiliary distros: 
	cout << "Finding auxiliary distros: " << endl; 
	fInName.str(""); 
	fInName << dataPathIn.str() << "nodesID_" << netName << ".csv"; 
	gNet.findDistros(fInName.str()); 
	cout << "Auxiliary distros found! " << endl; 

	// Writing auxiliary info to files in case the code crashes or we need to interrupt it: 
	cout << "Writing auxiliary info to files in case the code crashes or we need to interrupt it: " << endl; 
	fOutName.str(""); 
	fOutName << dataPathOut.str() << "hDistro_ref1.csv"; 
	gNet.saveHDistro(fOutName.str()); 
	fOutName.str(""); 
	fOutName << dataPathOut.str() << "hDistro_ref2.csv"; 
	gNet.saveHDistro_(fOutName.str()); 
	fOutName.str(""); 
	fOutName << dataPathOut.str() << "distantDistro.csv"; 
	gNet.saveDistantDistro(fOutName.str()); 
	fOutName.str(""); 
	fOutName << dataPathOut.str() << "diffsToRef_" << refID << ".csv"; 
	gNet.writeDiffToRef1(fOutName.str()); 
	cout << "Auxiliary info written to files in case the code crashes or we need to interrupt it! " << endl; 

	// Reporting execution time: 
	auto stop1 = high_resolution_clock::now(); 
	auto duration = duration_cast<seconds>(stop1 - start); 
	cout << "Computing distro time: " << duration.count() << endl; 

	// Finding Grand Network from detailed distances to reference: 
	cout << "Finding Grand Network from detailed distances to reference: " << endl; 
	gNet.findConnections(); 
	gNet.sortGConnections(); 
	cout << "Grand Network found from detailed distances to reference! " << endl; 

	// Saving Grand Network: 
	cout << "Saving Grand Network: " << endl; 
	fOutName.str(""); 
	fOutName << dataPathOut.str() << "grandNetwork.csv"; 
	gNet.saveGNet(fOutName.str()); 
	fOutName.str(""); 
	fOutName << dataPathOut.str() << "edges.csv"; 
	gNet.saveEdges(fOutName.str()); 
	cout << "Grand Network saved! " << endl; 

	auto stop2 = high_resolution_clock::now(); 
	duration = duration_cast<seconds>(stop2 - stop1); 
	cout << "Computing connections time: " << duration.count() << endl; 

	return 0;
}



