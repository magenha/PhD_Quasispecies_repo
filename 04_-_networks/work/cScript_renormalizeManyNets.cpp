/*

	cScript_renormalizeManyNets.cpp: 

		This script generates several renormalized networks with some common theme. For example, networks
		missing 10 consecutive nucleotides, or missing 10 nucleotides in random positions. The script
		reads a gNet, then generates a copy of it that is uses to produce a renormalized version. The
		renormalized network is stored in a folder, the folder named accordingly, and the process starts
		all over. In each case, the script also stores which nucleotides are omitted and which are kept. 

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
	string netName="r1"; 
	ostringstream dataPathIn, netPathIn, dataPathOut, sysStr, fInName, fOutName;
	dataPathIn.str(""); 
	dataPathIn << "/home/brigan/Desktop/Research_CNB/Bio/DeepSequencing/Data/GrandNetwork/Reg1/"; 
	cout << "Variables defined and paths loaded! " << endl; 

	// Generating several renormalized networks: 
	srand(23497); 
	int nRenormalizations = 10, nErase = 250; 
	for (int iR=0; iR<nRenormalizations; iR++){

		// Start chrono: 
		auto start = high_resolution_clock::now(); 
		// Loading network and assistant variables: 
		gN::GNet rNet; 
		fInName.str(""); 
		fInName << dataPathIn.str() << "metadata_" << netName << ".csv"; 
		rNet.loadMetadata(fInName.str()); 
		netPathIn.str(""); 
		netPathIn << dataPathIn.str() << "GN/"; 
		cout << "Loading network and assistant variables: " << endl; 
		fInName.str(""); 
		fInName << netPathIn.str() << "grandNetwork.csv"; 
		rNet.loadGNet(fInName.str());
		fInName.str("") ; 
		fInName << netPathIn.str() << "hDistro_ref1.csv"; 
		rNet.loadHDistro(fInName.str()); 
		fInName.str(""); 
		fInName << netPathIn.str() << "diffsToRef_" << rNet.getRefID() << ".csv"; 
		rNet.readDiffToRef1(fInName.str()); 

		auto stop0 = high_resolution_clock::now(); 
		auto duration = duration_cast<seconds>(stop0 - start); 
		cout << "Network and assistant variables loaded in " << duration.count() << "seconds!" << endl; 

		// Creating out-paths: 
		cout << "Creating out-paths for network " << iR << ": " << endl; 
		dataPathOut.str(""); 
		dataPathOut << dataPathIn.str() << "RN_" << iR << "_nE=" << nErase << "/"; 
		sysStr.str(""); 
		sysStr << "rm -r " << dataPathOut.str(); 
		system(sysStr.str().c_str()); 
		sysStr.str(""); 
		sysStr << "mkdir " << dataPathOut.str(); 
		system(sysStr.str().c_str()); 
		cout << "Out-paths created! " << endl; 

		// Defining renormalization protocol: 
		cout << "Defining " << iR << "-th renormalization protocol: " << endl; 
		vector<int> toErase; 
		fOutName.str(""); 
		fOutName << dataPathOut.str() << "protocol"; 
		toErase = h::rProtocol_nConsecutive(nErase, rNet.getSeqLength(), fOutName.str()); 
		cout << "Renormalization protocol defined. " << endl; 

		cout << "toErase: " << endl << "\t"; 
		for (int i=0; i<toErase.size(); i++){
			cout << toErase[i] << ", "; 
		}
		cout << endl; 

		// Renormalizing distros: 
		cout << "Renormalizing distros of network " << iR << ": " << endl; 
		rNet.renormalizeDistros(toErase); 
		auto stop1 = high_resolution_clock::now(); 
		duration = duration_cast<seconds>(stop1 - stop0); 
		cout << "Distros renormalized in " << duration.count() << " seconds. " << endl; 

		// Computing node catalog: 
		cout << "Computing node catalog for network " << iR << ": " << endl; 
		rNet.computeNodeCatalog(); 
		fOutName.str(""); 
		fOutName << dataPathOut.str() << "nodeCatalog.csv"; 
		rNet.saveNodeCatalog(fOutName.str()); 
		fOutName.str(""); 
		fOutName << dataPathOut.str() << "invertedNodeCatalog.csv"; 
		rNet.saveInvertedNodeCatalog(fOutName.str()); 
		auto stop2 = high_resolution_clock::now(); 
		duration = duration_cast<seconds>(stop2 - stop1); 
		cout << "Node catalog computed in " << duration.count() << " seconds. " << endl; 

		// Computing distros in renormalized space: 
		cout << "Computing distros in renormalized space for network " << iR << ": " << endl; 
		rNet.computeRepresentativeDistros(); 
		auto stop3 = high_resolution_clock::now(); 
		duration = duration_cast<seconds>(stop3 - stop2); 
		cout << "Distros in renormalized space computed in " << duration.count() << "seconds. " << endl; 

		// Computing connections in renormalized network: 
		cout << "Computing connections for renormalized network " << iR << ": " << endl; 
		rNet.findReConnections(); 
		rNet.sortGConnections(); 
		auto stop4 = high_resolution_clock::now(); 
		duration = duration_cast<seconds>(stop4 - stop3); 
		cout << "Connections of renormalized network computed in " << duration.count() << endl; 

		// Saving Renormalized Network: 
		cout << "Saving Renormalized Network " << iR << ": " << endl; 
		fOutName.str(""); 
		fOutName << dataPathOut.str() << "renormalizedNetwork.csv"; 
		rNet.saveGNet(fOutName.str()); 
		fOutName.str(""); 
		fOutName << dataPathOut.str() << "edges.csv"; 
		rNet.saveEdges(fOutName.str()); 
		cout << "Renormalized Network saved! " << endl; 

		// Writing distros to file: 
		cout << "Writing distros to file for network " << iR << ": " << endl; 
		fOutName.str(""); 
		fOutName << dataPathOut.str() << "hDistro_ref1.csv"; 
		rNet.saveHDistro(fOutName.str()); 
		fOutName.str(""); 
		fOutName << dataPathOut.str() << "distantDistro.csv"; 
		rNet.saveDistantDistro(fOutName.str()); 
		fOutName.str(""); 
		fOutName << dataPathOut.str() << "diffsToRef_" << rNet.getRefID() << ".csv"; 
		rNet.writeDiffToRef1(fOutName.str()); 
		cout << "Distros written to file! " << endl; 

	}

	return 0;
}



