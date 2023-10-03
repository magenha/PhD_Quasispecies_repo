/* 

	gNets.h: 

		This file implements genotype networks. 

*/

// Imports: 
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib> 
#include <vector>
#include <random>
#include <algorithm>
#include "helper.h"
using namespace std; 
// std::default_random_engine kMapGenerator; 


namespace gN{

	/*

		gNet namespace. 

	*/ 

	class GNet{
		/* GNet class: 

				This class defines objects of the kind of genotype networks. 

				Variables: 
					>> nGenotypes: Number of genotypes in the network. 
					>> nAllGenotypes: Number of genotypes in the original network. Only used after RG. 
					>> seqLength: Length of the sequences with which the current network works. 
					>> hThreshold: Threshold of Hamming distance below which two genotypes are considered linked. 
					>> hDistro: Distribution of Hamming distances to a reference sequence (usually the first one). 
					>> refID1, refID2: Indices of reference sequences. 
						> In some earlier versions, refID1=0 and refID2 points to the sequence more different to refID1. 
						> In a later version, refID1 should be the most common sequence across all files. 
							- This would minimize the amount of comparisons needed. 
						> This is used to pinpoint the Hamming distance to two references, thus increasing the constraints that must be
						satisfied before two sequences are tested in full length for connectivity. 
					>> hDistro_: Distro of Hamming distances to the auxiliary sequence. 
					>> gConnections: Array of vectors storing the genotype IDs to which each genotype is connected. 
					>> gID: Vector containing the IDs of each genotype in the network. 
						> This is needed for partial networks. 
						> We want to avoid loading all sequences in memory space for this code. 
						> Therefore, we will only load one or two sequences if needed. 
					>> distantSeqs, distantSeqs_: These arrays of vectors store the sequences that are a distance
					d^0 from ref 1: 		
						> This allows us to loop over distances, instead of over all sequences, when doing the
						comparison. 
					>> iDiffToRef1: Array of vectors containing the positions at which each sequence differs from reference 1. 
					>> diffToRef1: Nucleotides in the locations at which each sequence differs from reference 1. 
					>> whereMutations: Array of vectors storing the location along the sequence at which two connected genotypes differ. 
					>> nodeCatalog: This vector of vectors stores what sequences collapse into each other when
					some loci are ignored. This is only used for renormalized networks. 
					>> invertedNodeCatalog: This stores what renormalized node is each node mapped into. This way,
					it is very easy to retrieve such information. 

		*/ 

			int nGenotypes, nAllGenotypes, seqLength, hThreshold, *hDistro, refID1, refID2, *hDistro_; 
			vector<int> *gConnections, *whereMutations, gID, *distantSeqs, *distantSeqs_, *iDiffToRef1; 
			vector<char> *diffToRef1; 
			vector<vector<int>> nodeCatalog; 
			int *invertedNodeCatalog; 
		public: 
			// Init functions: 
			GNet(), GNet(string, int); 
			void findNGenotypes(string), findConnections(string); 
			void findHammingDistro(string, int), findHammingDistro_(string); 
			void findConnectionsFrom1Ref(string), findConnectionsFrom2Refs(string); 
			void findDistantDistro(string), findDistantDistro_(string), findConnectionsFromDistants(string); 
			void findDistros(string), findConnections(), findConectionsAndWhereMutations(); 
			void findReConnections(); 

			// Set functions: 
			void setNGenotypes(int), setHThreshold(int), setGConnections(vector<int> *&), setGID(vector<int> &); 
			void setRefID(int, int), setSeqLength(int), setSeqLengthFromFile(string); 

			// Get functions: 
			int getNGenotypes(), getNAllGenotypes(), getSeqLength(), getHThreshold(); 
			int getGConnection(int, int), getGID(int), getRefID(int); 
			vector<int> &getGConnections(int), *&getGConnections(), &getGID(); 
			vector<int> &getWhereMutations(int), *&getWhereMutations(); 
			int getHDistro(int), *&getHDistro(), getHDistro_(int), *&getHDistro_(); 
			vector<int> *&getDistantDistro(), &getDistantDistro(int), *&getDistantDistro_(), &getDistantDistro_(int); 
			vector<int> *&getIDiffToRef1(), &getIDiffToRef1(int); 
			vector<char> *&getDiffToRef1(), &getDiffToRef1(int); 
			vector<vector<int>> &getNodeCatalog(); 
			vector<int> &getNodeCatalog(int); 
			int *&getInvertedNodeCatalog(), getInvertedNodeCatalog(int); 

			// I/O functions: 
			void loadMetadata(string, int), readMetadata(string, int); 
			void saveGNet(string), writeGNet(string), loadGNet(string), readGNet(string); 
			void saveEdges(string), writeEdges(string), loadEdges(string), readEdges(string); 
			void saveHDistro(string), writeHDistro(string), loadHDistro(string), readHDistro(string); 
			void saveHDistro_(string), writeHDistro_(string), loadHDistro_(string), readHDistro_(string); 
			void saveDistantDistro(string), writeDistantDistro(string), loadDistantDistro(string), readDistantDistro(string); 
			void saveDiffToRef1(string), writeDiffToRef1(string), loadDiffToRef1(string), readDiffToRef1(string); 
			void saveWhereMutations(string), writeWhereMutations(string), loadWhereMutations(string), readWhereMutations(string); 
			void saveNodeCatalog(string), writeNodeCatalog(string), loadNodeCatalog(string), readNodeCatalog(string); 
			void saveInvertedNodeCatalog(string), writeInvertedNodeCatalog(string), loadInvertedNodeCatalog(string), readInvertedNodeCatalog(string); 

			// Functional functions: 
			void sortGConnections(); 
			void renormalizeDistros(vector<int> &), computeNodeCatalog(), computeRepresentativeDistros(); 

	}; 



	///////////////////////////////
	// 
	//  GNet functions: 

		///////////////////////////////
		// 
		//  Init Functions: 

		GNet::GNet(){return;} 

		GNet::GNet(string fInName, int hThreshold=1){
			/* Init function for objects of GNet class: 

				This function builds an object of class GNet. Therefore: 
					i) Reads a sorted list of unique genotypes from fIn. 
						> Genotypes are sorted, so the ID is the order in which they are read. 
					ii) Cheks, pair by pair, whether genotypes are below the threshold Hamming distance. 
					iii) Builds object gConnections that stores what IGs each genotype is connected to. 

				Inputs: 
					>> fInName: File from which genotypes can be read. 
					>> hThreshold: Maximum Hamming distance that still considers two genotypes linked. 

			*/ 

			findNGenotypes(fInName); 
			findConnections(fInName); 

			return; 
		}

		void GNet::findNGenotypes(string fInName){
			/* findNGenotypes function: 

				This function finds out the number of genotypes in the network. Therefore, it counts the number of lines in fIn. 

				Inputs: 
					>> fIn: File in which sorted genotypes have been stored, and from which the number of genotypes will be read. 

			*/ 

			// Word count to dummy file: 
			ostringstream strOut;
			strOut.str(""); 
			strOut << "wc -l < " << fInName.c_str() << " > temp.csv"; 
			system(strOut.str().c_str()); 

			// Dummy file content to integer: 
			strOut.str(""); 
			strOut << "temp.csv"; 
			ifstream fIn(strOut.str().c_str()); 
			string line; 
			getline(fIn, line); 
			setNGenotypes(h::getIntFromString(line)); 

			return; 
		}

		void GNet::findConnections(string fInName){
			/* findConnections function: 

				This function builds the genotype network from a given file. Therefore: 
					i) It loads one by one all couples of sequences possible. 
					ii) If runs the h::areConnected() function. 
					iii) If connected, the corresponding connections are added to the gConnections object. 
						> Note that this object has to be properly initialized. 

				ACHUNG!! This function is too slow for very large networks. We need an alternative. 

				Inputs: 
					>> fInName: File from which the network must be read. 

			*/ 

			// Initialize structure to store connections. 
			gConnections = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); // This is absurdly necessary... 
			}

			// Before starting to compare: 
			string seq1, seq2; 
			seq1 = h::loadGenotype(0, fInName); 
			seqLength = seq1.size(); 	// This way we don't need to compute the length repeatedly! 

			// Loop over all sequences: 
			// 	ACH! Last sequence is not visited in the main loop because it has already been visited by all others. 
			for (int iSeq1=0; iSeq1<nGenotypes-1; iSeq1++){
				cout << "Finding for: " << iSeq1 << endl; 
				seq1 = h::loadGenotype(iSeq1, fInName); 	// ACHTUNG!! Note that this needs to be adjusted if IDs change! 

				// Loop over sequences ahead of seq1: 
				for (int iSeq2=iSeq1+1; iSeq2<nGenotypes; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					if (h::areConnected(seqLength, seq1, seq2)){
						gConnections[iSeq1].push_back(iSeq2); 
						gConnections[iSeq2].push_back(iSeq1); 
					}
				}
			}

			return; 
		}

		// Functions to find network from references based on hDistros: 

		void GNet::findHammingDistro(string fInName, int iRefSeq=0){
			/* findHammingDistro function: 

				This function computes the Hamming distances of each sequence to a reference sequence (usually the first one in the
				provided file, fInName). 

				Inputs: 
					>> fInName: File from which sequences can be read. 
					>> iRefSeq=0: Reference sequence, usually the first one in the file. 

			*/ 

			// Initialize structure to store the distances: 
			hDistro = new int[nGenotypes]; 
			hDistro[iRefSeq]=0; 

			// Before starting to compare: 
			string seq1, seq2; 
			seq1 = h::loadGenotype(iRefSeq, fInName); 
			seqLength = seq1.size(); 	// This way we don't need to compute the length repeatedly! 

			// Loop over all other sequences: 
			// 	ACHTUNG!! By separating these two possibilities we spare ourselves a lot of IFs. 
			if (iRefSeq==0){
				for (int iSeq2=1; iSeq2<nGenotypes; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro[iSeq2] = h::hammingDistance(seqLength, seq1, seq2); 
				}
			}
			else{
				for (int iSeq2=0; iSeq2<nGenotypes; iSeq2++){
					if (iSeq2 != iRefSeq){
						seq2 = h::loadGenotype(iSeq2, fInName); 
						hDistro[iSeq2] = h::hammingDistance(seqLength, seq1, seq2); 
					}
				}
			}

			return; 
		}

		void GNet::findHammingDistro_(string fInName){
			/* findHammingDistro function: 

				This function computes the distributions of Hamming distances to the first sequence and to a second reference
				sequence that is the one farther apart from the first reference. 

				Inputs: 
					>> fInName: File from which sequences can be read. 

			*/ 


			// Initialize auxiliary sequence: 
			int maxH = -1; 
			refID2 = -1; 

			// Initialize structure to store the distro: 
			hDistro = new int[nGenotypes]; 
			hDistro[refID1] = 0; 

			// Before starting to compare: 
			string seq1, seq2, seq1_; 
			seq1 = h::loadGenotype(refID1, fInName); 
			seqLength = seq1.size(); 	// This way we don't need to compute the length repeatedly! 

			// Building distro for first reference: 
			// 	We do not need to compare refID1 with itself. 
			//  To minimize the number of operations we split the comparison in two batches: 
			for (int iSeq2=0; iSeq2<refID1; iSeq2++){
				seq2 = h::loadGenotype(iSeq2, fInName); 
				hDistro[iSeq2] = h::hammingDistance(seqLength, seq1, seq2); 
				if (hDistro[iSeq2] > maxH){
					maxH = hDistro[iSeq2]; 
					refID2 = iSeq2; 
				}
			}
			for (int iSeq2=refID1+1; iSeq2<nGenotypes; iSeq2++){
				seq2 = h::loadGenotype(iSeq2, fInName); 
				hDistro[iSeq2] = h::hammingDistance(seqLength, seq1, seq2); 
				if (hDistro[iSeq2] > maxH){
					maxH = hDistro[iSeq2]; 
					refID2 = iSeq2; 
				}
			}


			// Loading second reference: 
			seq1_ = h::loadGenotype(refID2, fInName); 

			// Initialize structure to store the second distro: 
			hDistro_ = new int[nGenotypes]; 
			hDistro_[refID1] = hDistro[refID2]; 
			hDistro_[refID2] = 0; 

			// Splitting this in three steps as well to minimize number of operations: 
			if (refID1 < refID2){
				for (int iSeq2=0; iSeq2<refID1; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
				for (int iSeq2=refID1+1; iSeq2<refID2; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
				for (int iSeq2=refID2+1; iSeq2<nGenotypes; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
			}
			else{
				for (int iSeq2=0; iSeq2<refID2; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
				for (int iSeq2=refID2+1; iSeq2<refID1; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
				for (int iSeq2=refID1+1; iSeq2<nGenotypes; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
			}

			return; 
		}

		void GNet::findConnectionsFrom1Ref(string fInName){
			/* findConnectionsFrom1Ref function: 

				This function uses the Hamming distance to a reference sequence to guide the seach for connections in the genotype
				network. It is assumed that the object hDistro has been computed already. 

				Let d^0_i and d^0_j be the distances to the reference sequence of sequences i and j respectively. Then, the
				following possibilities can happen:  
					>> If ||d^0_i - d^0_j|| > 1: 
						> There might be min(d^0_i, d^0_j) common mutations; but there are still a number >=2 outstanding mutations, hence
						i and j cannot be connected. 
					>> If ||d^0_i - d^0_j|| = 1: 
						> For i and j to be connected, there must be min(d^0_i, d^0_j) common mutations. 
					>> If ||d^0_i - d^0_j|| = 0: 
						> All mutations must be in the same positions -- otherwise the distance between i and j will be at least 2. 
						> Among these, at most one position can differ. 

				The current version does not exploit all the information that can be obtained from the reference sequence. At the
				moment, we just see if ||d^0_i - d^0_j||<=1. If so, we proceed with a crude comparison between sequences i and j. 

				Inputs: 
					>> fInName: From which to read the sequences. 

			*/ 

			// Initialize structure to store connections. 
			gConnections = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); // This is absurdly necessary... 
			}

			// Before starting to compare: 
			string seq1, seq2; 
			seq1 = h::loadGenotype(0, fInName); 
			int seqLength; 
			seqLength = seq1.size(); 	// This way we don't need to compute the length repeatedly! 

			// First build connections for the reference sequence: 
			for (int iSeq2=1; iSeq2<nGenotypes; iSeq2++){
				if (hDistro[iSeq2]==1){
					gConnections[0].push_back(iSeq2); 
					gConnections[iSeq2].push_back(0); 					
				}
			} 

			// Second, compare sequence distances one by one: 
			for (int iSeq1=1; iSeq1<nGenotypes-1; iSeq1++){
				cout << "Comparing to: " << iSeq1 << endl; 
				seq1 = h::loadGenotype(iSeq1, fInName); 

				// Second loop over sequences: 
				for (int iSeq2=iSeq1+1; iSeq2<nGenotypes; iSeq2++){
					// If differences of distances to reference is smaller than 2, they might be connected: 
					// 	Note that this step can be improved upon if we store info about mutations w.r.t. reference sequence. 
					if (abs(hDistro[iSeq1] - hDistro[iSeq2]) < 2){
						// Then, load the second sequence: 
						seq2 = h::loadGenotype(iSeq2, fInName); 
						// And check if they are connected: 
						if (h::areConnected(seqLength, seq1, seq2)){
							gConnections[iSeq1].push_back(iSeq2); 
							gConnections[iSeq2].push_back(iSeq1); 
						}
					}
				}
			}

			return; 
		}

		void GNet::findConnectionsFrom2Refs(string fInName){
			/* findConnectionsFrom2Refs function: 

				This function uses the same strategy as the previous one, but now taking two reference points. This allows us to
				"triangulate" the position of genotypes in the network, which discards even more candidates and accelerates the
				process. 

				The current version does not exploit all the information that can be obtained from the reference sequence. At the
				moment, we just see if ||d^0_i - d^0_j||<=1. If so, we proceed with a crude comparison between sequences i and j. 

				Inputs: 
					>> fInName: From which to read the sequences. 

			*/ 

			// Initialize structure to store connections. 
			gConnections = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); // This is absurdly necessary... 
			}

			// Before starting to compare: 
			string seq1, seq2; 
			seq1 = h::loadGenotype(0, fInName); 
			int seqLength; 
			seqLength = seq1.size(); 	// This way we don't need to compute the length repeatedly! 

			// First build connections for the reference sequence: 
			for (int iSeq2=1; iSeq2<nGenotypes; iSeq2++){
				if (hDistro[iSeq2]==1){
					gConnections[0].push_back(iSeq2); 
					gConnections[iSeq2].push_back(0); 					
				}
			} 
			// We could build connections for the second reference sequence as well, but by now connections in the gConnections
			// object are sorted in growing order, so I've decided to keep this by now. We won't be gaining much anyway. 

			// Second, compare sequence distances one by one: 
			bool fLoaded; // Tells us whether seq1 has been loaded already: 
			for (int iSeq1=1; iSeq1<nGenotypes-1; iSeq1++){
				cout << "Comparing to: " << iSeq1 << endl; 

				// Second loop over sequences: 
				fLoaded = false; 
				for (int iSeq2=iSeq1+1; iSeq2<nGenotypes; iSeq2++){
					// If differences of distances to reference is smaller than 2, they might be connected: 
					// 	Note that this step can be improved upon if we store info about mutations w.r.t. reference sequence. 
					if ((abs(hDistro[iSeq1] - hDistro[iSeq2]) < 2) and (abs(hDistro_[iSeq1] - hDistro_[iSeq2]) < 2)){
						if (not fLoaded){
							seq1 = h::loadGenotype(iSeq1, fInName); 
							fLoaded = true; 
						}
						// Then, load the second sequence: 
						seq2 = h::loadGenotype(iSeq2, fInName); 
						// And check if they are connected: 
						if (h::areConnected(seqLength, seq1, seq2)){
							gConnections[iSeq1].push_back(iSeq2); 
							gConnections[iSeq2].push_back(iSeq1); 
						}
					}
				}
			}

			return; 
		}


		// Functions to find network from references based on distants: 

		void GNet::findDistantDistro(string fInName){
			/* findDistantDistro function: 

				This function finds the distro of Hamming distances and also keeps track of lists of sequences that lie a given
				distance away. This will shorten the network reconstruction, since this will allow us to loop over distances,
				instead of over sequences, to do the comparison. 

				This function only keeps track of distant sequences w.r.t. the first reference sequence. 

				Inputs: 
					>> fInName: File from which sequences can be read. 

			*/ 


			// Initialize auxiliary sequence: 
			int maxH = -1, dummyH; 
			refID2 = -1; 

			// Initialize structure to store the distro: 
			hDistro = new int[nGenotypes]; 
			hDistro[refID1] = 0; 

			// Before starting to compare: 
			string seq1, seq2, seq1_; 
			seq1 = h::loadGenotype(refID1, fInName); 
			seqLength = seq1.size(); 	// This way we don't need to compute the length repeatedly! 

			// Initialize structure to get the distant distros: 
			// 	Note there is a similar distro w.r.t. the second reference sequence, but this is not used at the moment. 
			distantSeqs = new vector<int>[seqLength+1]; 
			for (int iDistance=0; iDistance<seqLength+1; iDistance++){
				distantSeqs[iDistance].clear(); 
			}


			// Building distro for first reference: 
			// 	As before, we split this process in two to minimize the number of operations: 
			distantSeqs[0].push_back(refID1); // Adding itself. 
			for (int iSeq2=0; iSeq2<refID1; iSeq2++){
				seq2 = h::loadGenotype(iSeq2, fInName); 
				dummyH = h::hammingDistance(seqLength, seq1, seq2); 
				hDistro[iSeq2] = dummyH; 
				distantSeqs[dummyH].push_back(iSeq2); 
				if (dummyH > maxH){
					maxH = hDistro[iSeq2]; 
					refID2 = iSeq2; 
				}
			}
			for (int iSeq2=refID1+1; iSeq2<nGenotypes; iSeq2++){
				seq2 = h::loadGenotype(iSeq2, fInName); 
				dummyH = h::hammingDistance(seqLength, seq1, seq2); 
				hDistro[iSeq2] = dummyH; 
				distantSeqs[dummyH].push_back(iSeq2); 
				if (dummyH > maxH){
					maxH = hDistro[iSeq2]; 
					refID2 = iSeq2; 
				}
			}


			// Loading recond reference: 
			seq1_ = h::loadGenotype(refID2, fInName); 

			// Initialize structure to store the second distro: 
			hDistro_ = new int[nGenotypes]; 
			hDistro_[refID1] = hDistro[refID2]; 
			hDistro_[refID2] = 0; // Distance of second reference to itself. 

			// As in the function before, we can decompose this loop in three steps: 
			if (refID1 < refID2){
				for (int iSeq2=0; iSeq2<refID1; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
				for (int iSeq2=refID1+1; iSeq2<refID2; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
				for (int iSeq2=refID2+1; iSeq2<nGenotypes; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
			}
			else{
				for (int iSeq2=0; iSeq2<refID2; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
				for (int iSeq2=refID2+1; iSeq2<refID1; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
				for (int iSeq2=refID1+1; iSeq2<nGenotypes; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					hDistro_[iSeq2] = h::hammingDistance(seqLength, seq1_, seq2); 
				}
			}

			return; 
		}

		void GNet::findDistantDistro_(string fInName){
			/* findDistantDistro_ function: 

				This function finds the distro of Hamming distances and also keeps track of lists of sequences that lie a given
				distance away. This will shorten the network reconstruction, since this will allow us to loop over distances,
				instead of over sequences, to do the comparison. 

				ACHTUNG!! 
					This functino hasn't been tested because it's not used. 

				Inputs: 
					>> fInName: File from which sequences can be read. 

			*/ 


			// Initialize auxiliary sequence: 
			int maxH = -1, dummyH; 
			refID2 = -1; 

			// Initialize structure to store the distro: 
			hDistro = new int[nGenotypes]; 
			hDistro[refID1] = 0; 

			// Before starting to compare: 
			string seq1, seq2, seq1_; 
			seq1 = h::loadGenotype(refID1, fInName); 
			seqLength = seq1.size(); 	// This way we don't need to compute the length repeatedly! 

			// Initialize structure to get the distant distros: 
			// 	Note there is a similar distro w.r.t. the second reference sequence, but this is not used at the moment. 
			distantSeqs = new vector<int>[seqLength+1]; 
			distantSeqs_ = new vector<int>[seqLength+1]; 
			for (int iDistance=0; iDistance<seqLength+1; iDistance++){
				distantSeqs[iDistance].clear(); 
				distantSeqs_[iDistance].clear(); 
			}


			// Building distro for first reference: 
			// 	As before, we split this process in two to minimize the number of operations: 
			distantSeqs[0].push_back(refID1); // Adding itself. 
			for (int iSeq2=0; iSeq2<refID1; iSeq2++){
				seq2 = h::loadGenotype(iSeq2, fInName); 
				dummyH = h::hammingDistance(seqLength, seq1, seq2); 
				hDistro[iSeq2] = dummyH; 
				distantSeqs[dummyH].push_back(iSeq2); 
				if (dummyH > maxH){
					maxH = hDistro[iSeq2]; 
					refID2 = iSeq2; 
				}
			}
			for (int iSeq2=refID1+1; iSeq2<nGenotypes; iSeq2++){
				seq2 = h::loadGenotype(iSeq2, fInName); 
				dummyH = h::hammingDistance(seqLength, seq1, seq2); 
				hDistro[iSeq2] = dummyH; 
				distantSeqs[dummyH].push_back(iSeq2); 
				if (dummyH > maxH){
					maxH = hDistro[iSeq2]; 
					refID2 = iSeq2; 
				}
			}

			// Loading recond reference: 
			seq1_ = h::loadGenotype(refID2, fInName); 

			// Initialize structure to store the second distro: 
			hDistro_ = new int[nGenotypes]; 
			hDistro_[refID1] = hDistro[refID2]; 
			hDistro_[refID2] = 0; // Distance of second reference to itself. 

			// As in the function before, we can decompose this loop in three steps: 
			if (refID1 < refID2){
				for (int iSeq2=0; iSeq2<refID1; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					dummyH = h::hammingDistance(seqLength, seq1_, seq2); 
					hDistro_[iSeq2] = dummyH; 
					distantSeqs_[dummyH].push_back(iSeq2); 
				}
				for (int iSeq2=refID1+1; iSeq2<refID2; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					dummyH = h::hammingDistance(seqLength, seq1_, seq2); 
					hDistro_[iSeq2] = dummyH; 
					distantSeqs_[dummyH].push_back(iSeq2); 
				}
				for (int iSeq2=refID2+1; iSeq2<nGenotypes; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					dummyH = h::hammingDistance(seqLength, seq1_, seq2); 
					hDistro_[iSeq2] = dummyH; 
					distantSeqs_[dummyH].push_back(iSeq2); 
				}
			}
			else{
				for (int iSeq2=0; iSeq2<refID2; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					dummyH = h::hammingDistance(seqLength, seq1_, seq2); 
					hDistro_[iSeq2] = dummyH; 
					distantSeqs_[dummyH].push_back(iSeq2); 
				}
				for (int iSeq2=refID2+1; iSeq2<refID1; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					dummyH = h::hammingDistance(seqLength, seq1_, seq2); 
					hDistro_[iSeq2] = dummyH; 
					distantSeqs_[dummyH].push_back(iSeq2); 
				}
				for (int iSeq2=refID1+1; iSeq2<nGenotypes; iSeq2++){
					seq2 = h::loadGenotype(iSeq2, fInName); 
					dummyH = h::hammingDistance(seqLength, seq1_, seq2); 
					hDistro_[iSeq2] = dummyH; 
					distantSeqs_[dummyH].push_back(iSeq2); 
				}
			}

			return; 
		}

		void GNet::findConnectionsFromDistants(string fInName){
			/* findConnectionsFromDistants function: 

				Previous functions to find the connections looped over all sequences and made comparisons (either letter by letter
				or by comparing the distances of each candidate couple to the references). This function, instead, loops over
				distances to the reference sequences and only compares sequences that lie within the allowed distance bracket for
				sequences to be connected. 

				Inputs: 
					>> fInName: From which to read the sequences. 

			*/ 

			// Initialize structure to store connections. 
			gConnections = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); // This is absurdly necessary... 
			}

			// Before starting to compare: 
			string seq1, seq2; 
			seq1 = h::loadGenotype(0, fInName); 
			int seqLength; 
			seqLength = seq1.size(); 	// This way we don't need to compute the length repeatedly! 
										// This is also the length over which to loop over. 

			// First build connections for the reference sequence: 
			for (int iDSeq=0; iDSeq<distantSeqs[1].size(); iDSeq++){
				gConnections[0].push_back(distantSeqs[1][iDSeq]); 
				gConnections[distantSeqs[1][iDSeq]].push_back(0); 
			}
			// We could build connections for the second reference sequence as well, but by now connections in the gConnections
			// object are sorted in growing order, so I've decided to keep this by now. We won't be gaining much anyway. 


			// Looping over all possible distances: 
			vector<int> dummyDistantSeqs1, dummyDistantSeqs2; 
			dummyDistantSeqs1.clear(); 
			dummyDistantSeqs2.clear(); 
			int iSeq1, iSeq2; 
			for (int iDistance=1; iDistance<seqLength; iDistance++){
				cout << "Comparing to distance: " << iDistance << endl; 
				// Load sequence IDs that are iDistance or iDistance+1 far apart from each other -- these might be connected: 
				dummyDistantSeqs1 = distantSeqs[iDistance]; 
				dummyDistantSeqs2 = distantSeqs[iDistance+1]; 

				// Elements in dummyDistantSeqs2 will be recalled again as iDistance grows in one unit. 
				// To avoid cross-comparing them twice, we know cross compare elements within dummyDistantSeqs1 and elements from
				// dummyDistantSeqs1 with elements of dummyDistantSeqs2, but we do not cross-compare elements within
				// dummyDistantSeqs2. 
				for (int iDSeq1=0; iDSeq1<dummyDistantSeqs1.size(); iDSeq1++){
					// Now it is very likely that we will need to do at least one comparison, so load seq1 right away: 
					iSeq1 = dummyDistantSeqs1[iDSeq1]; 
					seq1 = h::loadGenotype(iSeq1, fInName); 

					// Second loop over equally-distant sequences: 
					for (int iDSeq2=iDSeq1+1; iDSeq2<dummyDistantSeqs1.size(); iDSeq2++){
						iSeq2 = dummyDistantSeqs1[iDSeq2]; 
						seq2 = h::loadGenotype(iSeq2, fInName); 

						// Can we discard them with info from second reference sequence? 
						if (abs(hDistro_[iSeq1]-hDistro_[iSeq2])<2){
							// Are connected? 
							if (h::areConnected(seqLength, seq1, seq2)){
								gConnections[iSeq1].push_back(iSeq2); 
								gConnections[iSeq2].push_back(iSeq1); 
							} 
						}
					}

					// Loop over elements with 1-step higher relative distance: 
					for (int iDSeq2=0; iDSeq2<dummyDistantSeqs2.size(); iDSeq2++){
						iSeq2 = dummyDistantSeqs2[iDSeq2]; 
						seq2 = h::loadGenotype(iSeq2, fInName); 

						// Can we discard them with info from second reference sequence? 
						if (abs(hDistro_[iSeq1] - hDistro_[iSeq2])<2){
							// Are connected? 
							if (h::areConnected(seqLength, seq1, seq2)){
								gConnections[iSeq1].push_back(iSeq2); 
								gConnections[iSeq2].push_back(iSeq1); 
							} 
						}
					}
				}
			}

			return; 
		}


		// Functions to find connections based on detailed distances to refID1: 

		void GNet::findDistros(string fInName){
			/* findDistros function: 

				This function finds all the elements needed to make the faster network reconstruction. This is,
				with respecto to refID1, this function: i) Finds out hDistro, ii) finds out the distro of
				distants, iii) finds out and stores the exact differences of each sequences. 

				Note that, again, the sequence furthest away from refID1 is computed within this function -- so
				it will be ignored if provided. 

				Inputs: 
					>> fInName: From which to read the sequences. 					

			*/


			// Initialize auxiliary sequence: 
			int maxH = -1, dummyH; 
			refID2 = -1; 
			vector<int> dummyPositions; 

			// Initialize structure to store the distro: 
			hDistro = new int[nGenotypes]; 

			// Before starting to compare: 
			string seq1, seq2, seq1_; 
			seq1 = h::loadGenotype(refID1, fInName); 

			// Initialize structure to get the distant distros: 
			// 	Note there is a similar distro w.r.t. the second reference sequence, but this is not used at the moment. 
			distantSeqs = new vector<int>[seqLength+1]; 
			for (int iDistance=0; iDistance<seqLength+1; iDistance++){
				distantSeqs[iDistance].clear(); 
			}

			// Initialize structures where the differences are stored: 
			iDiffToRef1 = new vector<int>[nGenotypes]; 
			diffToRef1 = new vector<char>[nGenotypes]; 
			for (int iSeq=0; iSeq<nGenotypes; iSeq++){
				iDiffToRef1[iSeq].clear(); 
				diffToRef1[iSeq].clear(); 
			}

			// Storing information that we already know (refID1 w.r.t. itself): 
			hDistro[refID1] = 0; 				// Stating that its distance to itself is 0. 
			distantSeqs[0].push_back(refID1); 	// Adding itself at 0 distance. 
												// iDiffToRef[refID1] and diffToRef1[refID1] are just empty vectors. 

			// Building distro for first reference: 
			// 	As before, we split this process in two to minimize the number of operations: 
			for (int iSeq2=0; iSeq2<refID1; iSeq2++){
				seq2 = h::loadGenotype(iSeq2, fInName); 
				dummyPositions = h::hammingDistanceAndDifferences(seqLength, seq1, seq2); 
				dummyH = dummyPositions[dummyPositions.size()-1]; 
				dummyPositions.pop_back(); 
				hDistro[iSeq2] = dummyH; 
				distantSeqs[dummyH].push_back(iSeq2); 
				if (dummyH > maxH){
					maxH = hDistro[iSeq2]; 
					refID2 = iSeq2; 
				}
				for (int i=0; i<dummyPositions.size(); i++){
					iDiffToRef1[iSeq2].push_back(dummyPositions[i]); 
					diffToRef1[iSeq2].push_back(seq2[dummyPositions[i]]); 
				}
			}
			for (int iSeq2=refID1+1; iSeq2<nGenotypes; iSeq2++){
				seq2 = h::loadGenotype(iSeq2, fInName); 
				dummyPositions = h::hammingDistanceAndDifferences(seqLength, seq1, seq2); 
				dummyH = dummyPositions[dummyPositions.size()-1]; 
				dummyPositions.pop_back(); 
				hDistro[iSeq2] = dummyH; 
				distantSeqs[dummyH].push_back(iSeq2); 
				if (dummyH > maxH){
					maxH = hDistro[iSeq2]; 
					refID2 = iSeq2; 
				}
				for (int i=0; i<dummyPositions.size(); i++){
					iDiffToRef1[iSeq2].push_back(dummyPositions[i]); 
					diffToRef1[iSeq2].push_back(seq2[dummyPositions[i]]); 
				}
			}


			// Loading recond reference: 
			seq1_ = h::loadGenotype(refID2, fInName); 

			// Initialize structure to store the second distro: 
			hDistro_ = new int[nGenotypes]; 
			hDistro_[refID1] = hDistro[refID2]; 
			hDistro_[refID2] = 0; // Distance of second reference to itself. 

			// As in the function before, we can decompose this loop in three steps: 
			if (refID1 < refID2){
				for (int iSeq2=0; iSeq2<refID1; iSeq2++){
					hDistro_[iSeq2] = h::hammingDistanceFromDifferences(iDiffToRef1[refID2], iDiffToRef1[iSeq2], diffToRef1[refID2], diffToRef1[iSeq2]); 
				}
				for (int iSeq2=refID1+1; iSeq2<refID2; iSeq2++){
					hDistro_[iSeq2] = h::hammingDistanceFromDifferences(iDiffToRef1[refID2], iDiffToRef1[iSeq2], diffToRef1[refID2], diffToRef1[iSeq2]); 
				}
				for (int iSeq2=refID2+1; iSeq2<nGenotypes; iSeq2++){
					hDistro_[iSeq2] = h::hammingDistanceFromDifferences(iDiffToRef1[refID2], iDiffToRef1[iSeq2], diffToRef1[refID2], diffToRef1[iSeq2]); 
				}
			}
			else{
				for (int iSeq2=0; iSeq2<refID2; iSeq2++){
					hDistro_[iSeq2] = h::hammingDistanceFromDifferences(iDiffToRef1[refID2], iDiffToRef1[iSeq2], diffToRef1[refID2], diffToRef1[iSeq2]); 
				}
				for (int iSeq2=refID2+1; iSeq2<refID1; iSeq2++){
					hDistro_[iSeq2] = h::hammingDistanceFromDifferences(iDiffToRef1[refID2], iDiffToRef1[iSeq2], diffToRef1[refID2], diffToRef1[iSeq2]); 
				}
				for (int iSeq2=refID1+1; iSeq2<nGenotypes; iSeq2++){
					hDistro_[iSeq2] = h::hammingDistanceFromDifferences(iDiffToRef1[refID2], iDiffToRef1[iSeq2], diffToRef1[refID2], diffToRef1[iSeq2]); 
				}
			}

			return; 
		}

		void GNet::findConnections(){
			/* findConnections function: 

				ACHTUNG!! 
					This function uses the same shortcuts as the next function (findConectionsAndWhereMutations()) to
					reconstruct the Grand Network. The only difference is that this function does not produce
					whereMutations[]. In the future, we should use the next function only, as whereMutations[] will
					be required. I keep this by now because I would like to debug the sorting of connections --
					which has not been done yet, but it is not indispensable. 

				This function uses three pieces of information to find out the genotype networks: 
					i) It only checks couples of sequences which distances satisfy the required relationships to
					refID1. 
					ii) It checks the distance of each of these sequences to refID2. 
					iii) It looks only at those positions where each sequence differes from refID1, leveraging info
					of both which nucleotides are different and whether they are different. 

			*/ 

			// Initialize structure to store connections. 
			gConnections = new vector<int>[nGenotypes]; 
			whereMutations = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); // This is absurdly necessary... 
				whereMutations[i].clear(); // Still, absurdly necessary... 
			}

			// First build connections for the reference sequence -- therefore, we loop over seqs at distance 1: 
			for (int iDSeq=0; iDSeq<distantSeqs[1].size(); iDSeq++){
				gConnections[refID1].push_back(distantSeqs[1][iDSeq]); // Therefore we access distantSeqs[1]. 
				gConnections[distantSeqs[1][iDSeq]].push_back(refID1); // To navigate indexes more easily. 
			}

			// Looping over all possible distances: 
			vector<int> dummyDistantSeqs1, dummyDistantSeqs2; 
			dummyDistantSeqs1.clear(); 
			dummyDistantSeqs2.clear(); 
			int iSeq1, iSeq2; 
			for (int iDistance=1; iDistance<seqLength; iDistance++){
				cout << "Comparing for distance: " << iDistance << endl; 

				// Load sequence IDs that are iDistance or iDistance+1 far apart from each other -- these might be connected: 
				dummyDistantSeqs1 = distantSeqs[iDistance]; 
				dummyDistantSeqs2 = distantSeqs[iDistance+1]; 

				// Elements in dummyDistantSeqs2 will be recalled again as iDistance grows in one unit.  To avoid cross-comparing
				// them twice, we first cross-compare all elements within dummyDistantSeqs1 and elements from dummyDistantSeqs1 with
				// elements of dummyDistantSeqs2, but we do not cross-compare elements within dummyDistantSeqs2. 
				for (int iDSeq1=0; iDSeq1<dummyDistantSeqs1.size(); iDSeq1++){
					// Our comparisons will now be based on iDiffToRef1 and diffToRef1, so we don't need to load whole seqs anymore. 
					iSeq1 = dummyDistantSeqs1[iDSeq1]; 

					// Second loop over equally-distant sequences: 
					for (int iDSeq2=iDSeq1+1; iDSeq2<dummyDistantSeqs1.size(); iDSeq2++){
						iSeq2 = dummyDistantSeqs1[iDSeq2]; 

						// Can we discard them with info from second reference sequence? 
						if (abs(hDistro_[iSeq1]-hDistro_[iSeq2])<2){
							// Are connected? It is worth reading the documentation for this function! 
							if (h::areConnected(iDiffToRef1[iSeq1], iDiffToRef1[iSeq2], diffToRef1[iSeq1], diffToRef1[iSeq2])){
								gConnections[iSeq1].push_back(iSeq2); 
								gConnections[iSeq2].push_back(iSeq1); 
							} 
						}
					}

					// Loop over elements with 1-step higher relative distance: 
					for (int iDSeq2=0; iDSeq2<dummyDistantSeqs2.size(); iDSeq2++){
						iSeq2 = dummyDistantSeqs2[iDSeq2]; 

						// Can we discard them with info from second reference sequence? 
						if (abs(hDistro_[iSeq1] - hDistro_[iSeq2])<2){
							// Are connected? It is worth reading the documentation for this function! 
							if (h::areConnected(iDiffToRef1[iSeq1], iDiffToRef1[iSeq2], diffToRef1[iSeq1], diffToRef1[iSeq2])){
								gConnections[iSeq1].push_back(iSeq2); 
								gConnections[iSeq2].push_back(iSeq1); 
							} 
						}
					}
				}
			}

			return; 
		}

		void GNet::findConectionsAndWhereMutations(){
			/* findConectionsAndWhereMutations function: 

				ACHTUGN!! 
					This function implements the same reconstruction of the GN as the previous one and, additionally, returns
					whereMutations[]. This function has been tested and works properly. 

				This function uses three pieces of information to find out the genotype networks: 
					i) It only checks couples of sequences which distances satisfy the required relationships to refID1. 
					ii) It checks the distance of each of these sequences to refID2. 
					iii) It looks only at those positions where each sequence differes from refID1, leveraging info of both which
					nucleotides are different and whether they are different. 

				Inputs: 
					>> fInName: From which to read the sequences. 

			*/ 

			// Initialize structure to store connections. 
			gConnections = new vector<int>[nGenotypes]; 
			whereMutations = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); // This is absurdly necessary... 
				whereMutations[i].clear(); // Still, absurdly necessary... 
			}

			// Variables to handle whereMutations: 
			int iMuPlus, thisWhereMutation; 
			bool bMuPlus; 

			// First build connections for the reference sequence -- therefore, we loop over seqs at distance 1: 
			for (int iDSeq=0; iDSeq<distantSeqs[1].size(); iDSeq++){
				gConnections[refID1].push_back(distantSeqs[1][iDSeq]); // Therefore we access distantSeqs[1]. 
				gConnections[distantSeqs[1][iDSeq]].push_back(refID1); // To navigate indexes more easily. 
				thisWhereMutation = iDiffToRef1[distantSeqs[1][iDSeq]][0]; 
				whereMutations[refID1].push_back(thisWhereMutation); 
				whereMutations[distantSeqs[1][iDSeq]].push_back(thisWhereMutation); 
			}

			// Looping over all possible distances: 
			vector<int> dummyDistantSeqs1, dummyDistantSeqs2; 
			dummyDistantSeqs1.clear(); 
			dummyDistantSeqs2.clear(); 
			int iSeq1, iSeq2; 
			for (int iDistance=1; iDistance<seqLength; iDistance++){
				cout << "Comparing for distance: " << iDistance << endl; 

				// Load sequence IDs that are iDistance or iDistance+1 far apart from each other -- these might be connected: 
				dummyDistantSeqs1 = distantSeqs[iDistance]; 
				dummyDistantSeqs2 = distantSeqs[iDistance+1]; 

				// Elements in dummyDistantSeqs2 will be recalled again as iDistance grows in one unit.  To avoid cross-comparing
				// them twice, we first cross-compare all elements within dummyDistantSeqs1 and elements from dummyDistantSeqs1 with
				// elements of dummyDistantSeqs2, but we do not cross-compare elements within dummyDistantSeqs2. 
				for (int iDSeq1=0; iDSeq1<dummyDistantSeqs1.size(); iDSeq1++){
					// Our comparisons will now be based on iDiffToRef1 and diffToRef1, so we don't need to load whole seqs anymore. 
					iSeq1 = dummyDistantSeqs1[iDSeq1]; 

					// Second loop over equally-distant sequences: 
					for (int iDSeq2=iDSeq1+1; iDSeq2<dummyDistantSeqs1.size(); iDSeq2++){
						iSeq2 = dummyDistantSeqs1[iDSeq2]; 

						// Can we discard them with info from second reference sequence? 
						if (abs(hDistro_[iSeq1]-hDistro_[iSeq2])<2){
							// Are connected? It is worth reading the documentation for this function! 
							iMuPlus = h::areConnected_(iDiffToRef1[iSeq1], iDiffToRef1[iSeq2], diffToRef1[iSeq1], diffToRef1[iSeq2]); 
							if (iMuPlus){
								gConnections[iSeq1].push_back(iSeq2); 
								gConnections[iSeq2].push_back(iSeq1); 
								whereMutations[iSeq1].push_back(iMuPlus-1); 
								whereMutations[iSeq2].push_back(iMuPlus-1); 
							} 
						}
					}

					// Loop over elements with 1-step higher relative distance: 
					for (int iDSeq2=0; iDSeq2<dummyDistantSeqs2.size(); iDSeq2++){
						iSeq2 = dummyDistantSeqs2[iDSeq2]; 

						// Can we discard them with info from second reference sequence? 
						if (abs(hDistro_[iSeq1] - hDistro_[iSeq2])<2){
							// Are connected? It is worth reading the documentation for this function! 
							iMuPlus = h::areConnected_(iDiffToRef1[iSeq1], iDiffToRef1[iSeq2], diffToRef1[iSeq1], diffToRef1[iSeq2]); 
							if (iMuPlus){
								gConnections[iSeq1].push_back(iSeq2); 
								gConnections[iSeq2].push_back(iSeq1); 
								whereMutations[iSeq1].push_back(iMuPlus-1); 
								whereMutations[iSeq2].push_back(iMuPlus-1); 
							} 
						}
					}
				}
			}

			return; 
		}

		void GNet::findReConnections(){
			/* findReConnections function: 

				ACHTUGN!!  

					This function implements the same reconstruction of the GN as the previous, but it does not
					make use of the distance to the second reference. This is thought of to reconstruct
					renormalized networks, in which the distance to the second reference has been lost. 

				This function uses three pieces of information to find out the genotype networks: 
					i) It only checks couples of sequences which distances satisfy the required relationships to refID1. 
					ii) It looks only at those positions where each sequence differes from refID1, leveraging info
					of both which nucleotides are different and whether they are different. 

			*/ 

			// Initialize structure to store connections. 
			gConnections = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); // This is absurdly necessary... 
			}

			// First build connections for the node that contains the reference sequence. 
			// Therefore, we loop over seqs at distance 1: 
			for (int iDSeq=0; iDSeq<distantSeqs[1].size(); iDSeq++){
				gConnections[refID1].push_back(distantSeqs[1][iDSeq]); // Therefore we access distantSeqs[1]. 
				gConnections[distantSeqs[1][iDSeq]].push_back(refID1); // To navigate indexes more easily. 
			}

			// Looping over all possible distances: 
			vector<int> dummyDistantSeqs1, dummyDistantSeqs2; 
			dummyDistantSeqs1.clear(); 
			dummyDistantSeqs2.clear(); 
			int iSeq1, iSeq2; 
			for (int iDistance=1; iDistance<seqLength; iDistance++){
				cout << "Comparing for distance: " << iDistance << endl; 

				// Load sequence IDs that are iDistance or iDistance+1 far apart from each other -- these might be connected: 
				dummyDistantSeqs1 = distantSeqs[iDistance]; 
				dummyDistantSeqs2 = distantSeqs[iDistance+1]; 

				// Elements in dummyDistantSeqs2 will be recalled again as iDistance grows in one unit.  To avoid cross-comparing
				// them twice, we first cross-compare all elements within dummyDistantSeqs1 and elements from dummyDistantSeqs1 with
				// elements of dummyDistantSeqs2, but we do not cross-compare elements within dummyDistantSeqs2. 
				for (int iDSeq1=0; iDSeq1<dummyDistantSeqs1.size(); iDSeq1++){
					// Our comparisons will now be based on iDiffToRef1 and diffToRef1, so we don't need to load whole seqs anymore. 
					iSeq1 = dummyDistantSeqs1[iDSeq1]; 

					// Second loop over equally-distant sequences: 
					for (int iDSeq2=iDSeq1+1; iDSeq2<dummyDistantSeqs1.size(); iDSeq2++){
						iSeq2 = dummyDistantSeqs1[iDSeq2]; 

						// Are connected? 
						if (h::areConnected(iDiffToRef1[iSeq1], iDiffToRef1[iSeq2], diffToRef1[iSeq1], diffToRef1[iSeq2])){
							gConnections[iSeq1].push_back(iSeq2); 
							gConnections[iSeq2].push_back(iSeq1); 
						} 
					}

					// Loop over elements with 1-step higher relative distance: 
					for (int iDSeq2=0; iDSeq2<dummyDistantSeqs2.size(); iDSeq2++){
						iSeq2 = dummyDistantSeqs2[iDSeq2]; 
						// Are connected? 
						if (h::areConnected(iDiffToRef1[iSeq1], iDiffToRef1[iSeq2], diffToRef1[iSeq1], diffToRef1[iSeq2])){
							gConnections[iSeq1].push_back(iSeq2); 
							gConnections[iSeq2].push_back(iSeq1); 
						} 
					}
				}
			}

			return; 
		}

		///////////////////////////////
		// 
		//  Set Functions: 

		void GNet::setNGenotypes(int nGenotypes_){
			/* setNGenotypes function: 

			*/ 

			nGenotypes = nGenotypes_; 

			return; 
		}

		void GNet::setHThreshold(int hThreshold_){
			/* setHThreshold function: 

			*/ 

			hThreshold = hThreshold_; 

			return; 
		}

		void GNet::setGConnections(vector<int> *&gConnections_){
			/* setGConnections function: 

			*/ 

			gConnections = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); 
				for (int j=0; j<gConnections_[i].size(); j++){
					gConnections[i].push_back(gConnections_[i][j]); 
				}
			}

			return; 
		}

		void GNet::setGID(vector<int> &gID_){
			/* setGID function: 

			*/ 

			gID.clear(); 
			for (int i=0; i<nGenotypes; i++){
				gID.push_back(gID_[i]); 
			}

			return; 
		}

		void GNet::setRefID(int refID1_=0, int refID2_=0){
			/* setRefID function: 

			*/ 

			refID1 = refID1_; 
			refID2 = refID2_; 

			return; 
		}

		void GNet::setSeqLength(int seqLength_){
			/* setSeqLength function: 

			*/ 

			seqLength = seqLength_; 

			return; 
		}

		void GNet::setSeqLengthFromFile(string fInName){
			/* setSeqLengthFromFile function: 

				This function reads the first line of the sequence file just to read the length of the sequences. 

				Inputs: 
					>> fInName: Name of the file where the sequences are stored. 

			*/

			ifstream fIn(fInName.c_str()); 
			string seq; 
			seq = h::loadGenotype(0, fInName); 
			seqLength = seq.size(); 

			return; 
		}


		///////////////////////////////
		// 
		//  Get Functions: }

		int GNet::getNGenotypes(){
			/* getNGenotypes function: 

				Returns: 
					<< nGenotypes: Number of genotypes in the genotype network. 

			*/ 

			return nGenotypes; 
		}

		int GNet::getNAllGenotypes(){
			/* getNAllGenotypes function: 

				Returns: 
					<< nAllGenotypes: Used only after renormalization. Stores the original # of genotypes. 

			*/

			return nAllGenotypes; 
		}

		int GNet::getSeqLength(){
			/* getSeqLength function: 

				Returns: 
					<< seqLength: Length of the sequence of this genotype network. 

			*/

			return seqLength; 
		}

		int GNet::getHThreshold(){
			/* getHThreshold function: 

				Returns: 
					<< hThreshold: Threshold below which we consider two genotypes connected. 

			*/ 

			return hThreshold; 
		}

		int GNet::getGConnection(int i, int j){
			/* getGConnection function: 

				Inputs: 
					>> i, j: Index to retrieve the desired edge. 

				Returns: 
					<< gConnections[i][j]: j-th connectino of the i-th genotype. 
						< ACHTUNG!! Make sure to call this within bonds. 

			*/ 

			return gConnections[i][j]; 
		}

		int GNet::getGID(int i){
			/* getGID function: 

				Inputs: 
					>> i: Index to retrieve ID from. 

				Returns: 
					<< gID[i]: ID of the i-th genotype in the network. 

			*/ 

			return gID[i]; 
		}

		int GNet::getRefID(int iID=0){
			/* getRefID fucntion: 

				Inputs: 
					>> iID=0: Indicates which reference ID to retrieve. 
						> false: Retrieves the first reference ID. 
						> true: Retrieves the second reference ID. 

				Returns: 
					<< refID1: If not(iID). 
					<< refID2: If iID. 

			*/

			if (not(iID)) return refID1; 

			return refID2; 
		}

		vector<int> &GNet::getGConnections(int i){
			/* getGConnections function: 

				Inputs: 
					>> i: Index to retrieve connectinos from. 

				Returns: 
					<< gConnections[i]: Vector with connections of the i-th genotype. 

			*/ 

			return gConnections[i]; 
		}

		vector<int> *&GNet::getGConnections(){
			/* getGConnections function: 

				Returns: 
					<< gConnections: Array of vectors with the connections. 

			*/ 

			return gConnections; 
		}

		vector<int> &GNet::getGID(){
			/* getGID function: 

				Returns: 
					<< gID: vector with the IDs of the genotypes in the network. 

			*/ 

			return gID; 
		}

		vector<int> &GNet::getWhereMutations(int i){
			/* getWhereMutations function: 

				Returns: 
					<< whereMutations[i]: Vector with the location of the mutations of all seqs. connected to seq. i. 

			*/ 

			return whereMutations[i]; 
		}

		vector<int> *&GNet::getWhereMutations(){
			/* getWhereMutations function: 

				Returns: 
					<< whereMutations: Array of vectors with the locations of the mutations of all seqs. connected to each seq. 

			*/ 

			return whereMutations; 
		}

		int *&GNet::getHDistro(){
			/* getHDistro funciton: 

				Returns: 
					<< hDistro: Distribution of Hammind distances to the first reference sequence. 

			*/ 

			return hDistro; 
		}

		int GNet::getHDistro(int iSeq){
			/* getHDistro function: 

				Inputs: 
					>> iSeq: Sequence whose Hamming distance to the first reference we wish to know. 

				Returns: 
					<< hDistro[iSeq]: Hamming distance of the iSeq-th sequence to the first reference sequence. 

			*/ 

			return hDistro[iSeq]; 
		}

		int *&GNet::getHDistro_(){
			/* getHDistro_ funciton: 

				Returns: 
					<< hDistro_: Distribution of Hammind distances to the second reference sequence. 

			*/ 

			return hDistro_; 
		}

		int GNet::getHDistro_(int iSeq){
			/* getHDistro_ function: 

				Inputs: 
					>> iSeq: Sequence whose Hamming distance to the second reference we wish to know. 

				Returns: 
					<< hDistro_[iSeq]: Hamming distance of the iSeq-th sequence to the second reference sequence. 

			*/ 

			return hDistro_[iSeq]; 
		}

		vector<int> *&GNet::getDistantDistro(){
			/* getDistantDistro function: 

				Returns: 
					<< distantSeqs: Array of vectors, each containing the IDs of sequences at a given distance of the first reference. 

			*/ 

			return distantSeqs; 
		}

		vector<int> &GNet::getDistantDistro(int iDistance){
			/* getDistantDistro function: 

				Inputs: 
					>> iDistance: Distance of which we wish to know what sequences lie this far apart from the first reference. 

				Returns: 
					<< distantSeqs[iDistance]: Vector of sequences that lie iDistance far apart from the first reference. 

			*/ 

			return distantSeqs[iDistance]; 
		}

		vector<int> *&GNet::getDistantDistro_(){
			/* getDistantDistro_ function: 

				Returns: 
					<< distantSeqs_: Array of vectors, each containing the IDs of sequences at a given distance of the second
					reference. 

			*/ 

			return distantSeqs_; 
		}

		vector<int> &GNet::getDistantDistro_(int iDistance){
			/* getDistantDistro_ function: 

				Inputs: 
					>> iDistance: Distance of which we wish to know what sequences lie this far apart from the second reference. 

				Returns: 
					<< distantSeqs_[iDistance]: Vector of sequences that lie iDistance far apart from the second reference. 

			*/ 

			return distantSeqs_[iDistance]; 
		}

		vector<int> *&GNet::getIDiffToRef1(){
			/* getIDiffToRef1 fuction: 

				Returns: 
					<< iDiffToRef1: Array of vectors containing the differenteces of each sequence to the first reference sequence. 

			*/ 

			return iDiffToRef1; 
		}

		vector<int> &GNet::getIDiffToRef1(int iSeq){
			/* getIDiffToRef1 fuction: 
				
				Inputs: 
					>> iSeq: Reference of which we wish to retrieve the differences to the reference sequence. 

				Returns: 
					<< iDiffToRef1[iSeq]: Vectors containing the differenteces of sequence iSeq to the first reference sequence. 

			*/ 

			return iDiffToRef1[iSeq]; 
		}

		vector<char> *&GNet::getDiffToRef1(){
			/* getDiffToRef1 function: 

				Returns: 
					<< diffToRef1: Array of vectors containing each sequence's nucleotides which are different to those of the
					reference sequence. 

			*/ 

			return diffToRef1; 
		}

		vector<char> &GNet::getDiffToRef1(int iSeq){
			/* getDiffToRef1 function: 

				Returns: 
					<< diffToRef1[iSeq]: Vector containing the nucleotides of sequence iSeq which are different to those of the
					reference sequence. 

			*/ 

			return diffToRef1[iSeq]; 
		}

		vector<vector<int>> &GNet::getNodeCatalog(){
			/* getNodeCatalog function: 

				Returns: 
					<< nodeCatalog: Vector of vectors summarizing which sequences are renormalized into each node. 

			*/ 

			return nodeCatalog; 
		}

		vector<int> &GNet::getNodeCatalog(int iNode){
			/* getNodeCatalog function: 

				Returns: 
					<< nodeCatalog[iNode]: Vector of integers summarizing which sequences are renormalized to node
					iNode. 

			*/ 

			return nodeCatalog[iNode]; 
		}

		int *&GNet::getInvertedNodeCatalog(){
			/* getInvertedCatalog function: 

				Returns: 
					<< invertedNodeCatalog: Array of integers that links each sequence in the original network into a
					node of a renormalized network. 

			*/

			return invertedNodeCatalog; 
		}

		int GNet::getInvertedNodeCatalog(int iSeq){
			/* getInvertedCatalog function: 
	
				Returns: 
					<< invertedNodeCatalog[iSeq]: Node to which the iSeq-th sequence has been linked through
					renormalization. 

			*/

			return invertedNodeCatalog[iSeq];
		}


		///////////////////////////////
		// 
		//  I/O Functions: 

		void GNet::loadMetadata(string fInName, int iVar=-1){
			/* loadMetadata function: 

				This functions loads relevant information for this network from a metadata file. This info might
				be superseeded later on. 

				Inputs: 
					>> fInName: Name of the file where the metadata will be read from. 
					>> iVar=-1: Indicates which variables need to be loaded. 
						> iVar=-1: All variables are loaded. 
						> iVar=1, 2, or 3: Only the first, second, or third variable are loaded. 

			*/


			ifstream fIn(fInName.c_str()); 
			string line; 

			if (iVar==-1){
				getline(fIn, line); 
				setNGenotypes(h::getIntFromString(line)); 
				getline(fIn, line); 
				setRefID(h::getIntFromString(line)); 
				getline(fIn, line); 
				setSeqLength(h::getIntFromString(line)); 
			}
			else{
				for (int i=0; i<iVar; i++){
					getline(fIn, line); 
				}
				if (iVar==1) setNGenotypes(h::getIntFromString(line)); 
				if (iVar==2) setRefID(h::getIntFromString(line)); 
				if (iVar==3) setSeqLength(h::getIntFromString(line)); 
			}

			return; 
		}


			void loadMetadata(string, int), readMetadata(string, int); 




		void GNet::saveGNet(string fOutName){
			/* saveGNet function: 

				ACHTUNG!! These functions will be deprecated, as the new GNet objects also store where the mutations took place!
				Use, instead, the same function ended with an underscore. 

				This function saves the genotype network as a list of numbers that allow us to recompose the network if loaded
				later. The format is as follows: 
					i) First, the number of genotypes in the network. 
					ii) For each genotype: 
						iia) The number of other genotypes it connects to. 
						iib) Each of the IDs of the connected genotypes. 

				We should add other implementations that store, e.g., a list of connections. 

				Inputs: 
					>> fOutName: Where the data will be saved. 

			*/ 

			ofstream fOut(fOutName.c_str()); 
			fOut << nGenotypes << endl; 
			for (int i=0; i<nGenotypes; i++){
				fOut << gConnections[i].size() << endl; 
				for (int j=0; j<gConnections[i].size(); j++){
					fOut << gConnections[i][j] << endl; 
				}
			}

			return; 
		}

		void GNet::writeGNet(string fOutName){
			/* writeGNet function: 

				ACHTUNG!! These functions will be deprecated, as the new GNet objects also store where the mutations took place!
				Use, instead, the same function ended with an underscore. 

				This function wraps the previous one so that I can call either indistinctly. 

			*/ 

			saveGNet(fOutName); 

			return; 
		}

		void GNet::loadGNet(string fInName){
			/* loadGNet function: 

				ACHTUNG!! These functions will be deprecated, as the new GNet objects also store where the mutations took place!
				Use, instead, the same function ended with an underscore. 

				This function loads a network that has been stored in the same format as that produced by saveGNet(). 

				Inputs: 
					>> fInName: From which to read the network. 

			*/ 

			ifstream fIn(fInName.c_str()); 
			string line; 
			int nConnections; 

			getline(fIn, line); 
			setNGenotypes(h::getIntFromString(line)); 
			cout << line << endl; 
			gConnections = new vector<int>[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				gConnections[i].clear(); 
				getline(fIn, line); 
				nConnections = h::getIntFromString(line); 
				for (int j=0; j<nConnections; j++){
					getline(fIn, line); 
					gConnections[i].push_back(h::getIntFromString(line)); 
				}
			}

			return; 
		}

		void GNet::readGNet(string fInName){
			/* readGNet function: 

				ACHTUNG!! These functions will be deprecated, as the new GNet objects also store where the mutations took place!
				Use, instead, the same function ended with an underscore. 

				To have both options. 

			*/ 

			loadGNet(fInName); 

			return; 
		}

		void GNet::saveEdges(string fOutName){
			/* saveEdges function: 

				This function saves the edges in the network as sets of pairs of nodes. This should contain all the information
				needed to reconstruct the network with python's networkx, except for isolated genotypes. These must be loaded from a
				list of nodes -- but they are arguably of little importance. 

				Inputs: 
					>> fOutName: Name of the file on which we wish to store the edges. 

			*/ 

			// Open file: 
			ofstream fOut(fOutName.c_str()); 

			// Loop over all nodes: 
			for (int i=0; i<nGenotypes; i++){
				// Loop over nodes connected to each node: 
				for (int j=0; j<gConnections[i].size(); j++){
					// To only print edges once: 
					if (gConnections[i][j] > i) fOut << i << ", " << gConnections[i][j] << endl; 
				}
			}
		}

		void GNet::writeEdges(string fOutName){
			/* writeEdges function: 

				To have both available. 

			*/ 

			saveEdges(fOutName); 

			return; 
		}

		void GNet::saveHDistro(string fOutName){
			/* saveHDistro function: 

				This function saves the distribution of Hamming distances. 

				Inputs: 
					>> fOutName: Name of the file on which we wish to store the edges. 

			*/ 


			// Open file: 
			ofstream fOut(fOutName.c_str()); 

			for (int i = 0; i<nGenotypes; i++){
				fOut << hDistro[i] << endl; 
			}

			return; 
		}

		void GNet::writeHDistro(string fOutName){
			/* writeHDistro function: 

				To have both available. 

			*/ 

			saveHDistro(fOutName); 

			return; 
		}

		void GNet::loadHDistro(string fInName){
			/* loadHDistro function: 

				This function loads the distro of distances to the first reference. 

				Inputs: 
					>> fInName: Where the distro has been previously stored. 

			*/ 

			// Variables to load data: 
			ifstream fIn(fInName.c_str()); 
			string line; 

			hDistro = new int[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				getline(fIn, line); 
				hDistro[i] = h::getIntFromString(line); 
			}

			return; 
		}

		void GNet::readHDistro(string fInName){
			/* readHDistro function: 

				To have both options available. 

			*/

			loadHDistro(fInName); 

			return; 
		}

		void GNet::saveHDistro_(string fOutName){
			/* saveHDistro_ function: 

				This function saves the distribution of Hamming distances. 

				Inputs: 
					>> fOutName: Name of the file on which we wish to store the edges. 

			*/ 


			// Open file: 
			ofstream fOut(fOutName.c_str()); 

			for (int i = 0; i<nGenotypes; i++){
				fOut << hDistro_[i] << endl; 
			}

			return; 
		}

		void GNet::writeHDistro_(string fOutName){
			/* writeHDistro_ function: 

				To have both available. 

			*/ 

			saveHDistro_(fOutName); 

			return; 
		}

		void GNet::loadHDistro_(string fInName){
			/* loadHDistro_ function: 

				This function loads the distro of distances to the second reference. 

				Inputs: 
					>> fInName: Where the distro has been previously stored. 

			*/ 

			// Variables to load data: 
			ifstream fIn(fInName.c_str()); 
			string line; 

			hDistro_ = new int[nGenotypes]; 
			for (int i=0; i<nGenotypes; i++){
				getline(fIn, line); 
				hDistro_[i] = h::getIntFromString(line); 
			}

			return; 
		}

		void GNet::readHDistro_(string fInName){
			/* readHDistro_ function: 

				To have both options available. 

			*/

			loadHDistro_(fInName); 

			return; 
		}

		void GNet::saveDistantDistro(string fOutName){
			/* saveDistantDistro function: 

				This function saves the distantSeqs objects that contains the relation of elements at each given distance. 

				Inputs: 
					>> fOutName: Name of the file on which we wish to store the object. 

			*/ 


			// Open file: 
			ofstream fOut(fOutName.c_str()); 

			fOut << seqLength << endl; 
			for (int iDist=0; iDist<seqLength+1; iDist++){
				fOut << distantSeqs[iDist].size() << endl; 
				for (int iDSeq=0; iDSeq<distantSeqs[iDist].size(); iDSeq++){
					fOut << distantSeqs[iDist][iDSeq] << endl; 
				}
			}

			return; 
		}

		void GNet::writeDistantDistro(string fOutName){
			/* writeDistantDistro function 

				To have both options available. 

			*/ 

			saveDistantDistro(fOutName); 

			return; 
		}

		void GNet::loadDistantDistro(string fInName){
			/* loadDistantDistro function: 

				This function loads the position of mutations (and what mutations they are) with respect to the first reference. 

				Inputs: 
					>> fInName: File from which we will read said differences. 

			*/

			// Variables to load data: 
			ifstream fIn(fInName.c_str()); 
			string line; 
			int thisNSeqs; 


			getline(fIn, line); 
			seqLength = h::getIntFromString(line); 
			distantSeqs = new vector<int>[seqLength+1]; 

			for (int iLen=0; iLen<seqLength+1; iLen++){
				distantSeqs[iLen].clear(); 
				getline(fIn, line); 
				thisNSeqs = h::getIntFromString(line); 
				for (int iSeq=0; iSeq<thisNSeqs; iSeq++){
					getline(fIn, line); 
					distantSeqs[iLen].push_back(h::getIntFromString(line)); 
				}
			}

			return; 
		}

		void GNet::readDistantDistro(string fInName){
			/* readDistantDistro function: 

				To have both options available. 

			*/

			loadDistantDistro(fInName); 

			return; 
		}

		void GNet::saveDiffToRef1(string fOutName){
			/* saveDiffToRef1 function: 

				This function saves the structures iDiffToRef1 and diffToRef1 in a file. 

				Inputs: 
					>> fOutName: Name of the file on which we wish to store the object. 

			*/ 


			// Open file: 
			ofstream fOut(fOutName.c_str()); 

			// Loop over all sequences: 
			for (int iSeq=0; iSeq<nGenotypes; iSeq++){
				fOut << iDiffToRef1[iSeq].size() << endl; 
				for (int iDiff=0; iDiff < iDiffToRef1[iSeq].size(); iDiff++){
					fOut << iDiffToRef1[iSeq][iDiff] << endl; 
					fOut << diffToRef1[iSeq][iDiff] << endl; 
				}
			}

			return; 
		}

		void GNet::writeDiffToRef1(string fOutName){
			/* writeDiffToRef1 function 

				To have both options available. 

			*/ 

			saveDiffToRef1(fOutName); 

			return; 
		}

		void GNet::loadDiffToRef1(string fInName){
			/* loadDiffToRef1 function: 

				This function loads the position and nucleotide content of all mutations w.r.t. ref 1. 

				Inputs: 
					>> fInName: File from which to read the data. 

			*/

			// Variables to load data: 
			ifstream fIn(fInName.c_str()); 
			string line; 
			int thisNMutations; 

			iDiffToRef1 = new vector<int>[nGenotypes]; 
			diffToRef1 = new vector<char>[nGenotypes]; 
			for (int iGenotype=0; iGenotype<nGenotypes; iGenotype++){
				iDiffToRef1[iGenotype].clear(); 
				diffToRef1[iGenotype].clear(); 
				getline(fIn, line); 
				thisNMutations = h::getIntFromString(line); 
				for (int iMutation=0; iMutation< thisNMutations; iMutation++){
					getline(fIn, line); 
					iDiffToRef1[iGenotype].push_back(h::getIntFromString(line)); 
					getline(fIn, line); 
					diffToRef1[iGenotype].push_back(line[0]); 
				}
			}

			return; 
		}

		void GNet::readDiffToRef1(string fInName){
			/* readDiffToRef1 function: 

				To have both options available. 

			*/

			loadDiffToRef1(fInName); 

			return; 
		}

		void GNet::saveWhereMutations(string fOutName){
			/* saveWhereMutations function: 

				This function saves the positions at which mutations happened. 

				Inputs: 
					>> fOutName: File where we will store the positions at which the mutations happened. 

			*/ 

			// Open file: 
			ofstream fOut(fOutName.c_str()); 

			// Loop over all sequences to store where mutations: 
			for (int iSeq=0; iSeq<nGenotypes; iSeq++){
				// Loop over size of whereMutations of this sequence: 
				for (int iWhere=0; iWhere<whereMutations[iSeq].size(); iWhere++){
					fOut << whereMutations[iSeq][iWhere] << endl; 
				}
			}

			return; 
		}

		void GNet::writeWhereMutations(string fOutName){
			/* writeWhereMutations function: 

				To have both options available. 

			*/ 

			saveWhereMutations(fOutName); 

			return; 
		}

		void GNet::loadWhereMutations(string fInName){
			/* loadWhereMutations function: 

				This function loads whereMutations for a given network. Note that this can only be called once the GNet object has
				been properly loaded, since it relies on its sizes, etc. Also, it is assumed that connections of each node are
				sorted in the same way in the gConnections and whereMutations objects. 

				Inputs: 
					>> fInName: File from which whereMutations will be read. 

			*/ 

			// Opening file and defining dummy variables: 
			ifstream fIn(fInName.c_str()); 
			string line; 
			whereMutations = new vector<int>[nGenotypes]; 

			// Loop over sequences: 
			for (int iSeq=0; iSeq<nGenotypes; iSeq++){
				whereMutations[iSeq].clear(); 	// Still absurdly necessary... 
				// Loop over number of connections for this sequence: 
				for (int iWhere=0; iWhere<gConnections[iSeq].size(); iWhere++){
					getline(fIn, line); 
					whereMutations[iSeq].push_back(h::getIntFromString(line)); 
				}
			}

			return; 
		}

		void GNet::readWhereMutations(string fInName){
			/* readWhereMutations function: 

				To have both options available. 

			*/

			loadWhereMutations(fInName); 

			return; 
		}

		void GNet::saveNodeCatalog(string fOutName){
			/* saveNodeCatalog function: 

				This function saves the catalog that tells which original sequences are contained within each
				renormalized node. 

				Inputs: 
					>> fOutName: File where we will store the nodeCatalog. 

			*/

			// Open file: 
			ofstream fOut(fOutName.c_str()); 

			fOut << nGenotypes << endl; 
			for (int iNode=0; iNode<nGenotypes; iNode++){
				fOut << nodeCatalog[iNode].size() << endl; 
				for (int iSeq=0; iSeq<nodeCatalog[iNode].size(); iSeq++){
					fOut << nodeCatalog[iNode][iSeq] << endl; 
				}
			}

			return; 
		}

		void GNet::writeNodeCatalog(string fOutName){
			/* writeNodeCatalog function: 

				To have both options available. 

			*/

			saveNodeCatalog(fOutName); 

			return; 
		}

		void GNet::loadNodeCatalog(string fInName){
			/* loadNodeCatalog function: 

				This function loads a catalog of nodes. 

				Inputs: 
					>> fInName: Name of the file where the catalog is stored. 

			*/

			// Opening file and defining dummy variables: 
			ifstream fIn(fInName.c_str()); 
			string line; 
			nodeCatalog.clear(); 
			int nSeqs, thisSeq; 

			getline(fIn, line); 
			nGenotypes = h::getIntFromString(line); 
			for (int iNode=0; iNode<nGenotypes; iNode++){
				nodeCatalog.push_back({}); 
				getline(fIn, line); 
				nSeqs = h::getIntFromString(line); 
				for (int iSeq=0; iSeq<nSeqs; iSeq++){
					getline(fIn, line); 
					thisSeq = h::getIntFromString(line); 
					nodeCatalog[iNode].push_back(thisSeq); 
				}
			}

			return; 
		}

		void GNet::readNodeCatalog(string fInName){
			/* readNodeCatalog function: 

				To have both options available. 

			*/

			loadNodeCatalog(fInName); 

			return; 
		}

		void GNet::saveInvertedNodeCatalog(string fOutName){
			/* saveInvertedNodeCatalog function: 

				This function saves the inverted catalog, which allows to quickly access what node each sequence
				has been mapped through renormalization. 

				Inputs: 
					>> fOutName: Where we will write the inverted catalog. 

			*/

			// Open file: 
			ofstream fOut(fOutName.c_str()); 

			fOut << nAllGenotypes << endl; 
			for (int iSeq=0; iSeq<nAllGenotypes; iSeq++){
				fOut << invertedNodeCatalog[iSeq] << endl; 
			}

			return; 
		}

		void GNet::writeInvertedNodeCatalog(string fOutName){
			/* writeInvertedNodeCatalog function: 

				To have both options available. 

			*/

			saveInvertedNodeCatalog(fOutName); 

			return; 
		}

		void GNet::loadInvertedNodeCatalog(string fInName){
			/* loadInvertedNodeCatalog function: 

				This function loads an inverted node catalog from a file. 

				Inputs: 
					>> fInName: From where the inverted node catalog will be read.

			*/

			// Opening file and defining dummy variables: 
			ifstream fIn(fInName.c_str()); 
			string line; 

			getline(fIn, line); 
			nAllGenotypes = h::getIntFromString(line); 
			invertedNodeCatalog = new int[nAllGenotypes]; 
			for (int iNode=0; iNode<nAllGenotypes; iNode++){
				getline(fIn, line); 
				invertedNodeCatalog[iNode] = h::getIntFromString(line); 
			}
			return; 
		}

		void GNet::readInvertedNodeCatalog(string fInName){
			/* readInvertedNodeCatalog function: 

				To have both options available. 

			*/

			loadInvertedNodeCatalog(fInName); 

			return; 
		}


		///////////////////////////////
		// 
		//  Functional Functions: 

		void GNet::sortGConnections(){
			/* sortGConnections function: 

				This function sorts the elements within each gConnections vector. 

			*/ 

			bool fSort=true; 
			int dummyConnection; 

			// Loop over all gConnections: 
			for (int iGenotype=1; iGenotype<nGenotypes; iGenotype++){
				sort(gConnections[iGenotype].begin(), gConnections[iGenotype].end()); 
			}

			return; 
		}

		void GNet::renormalizeDistros(vector<int> &toIgnore){
			/* renormalizeDistros function: 

				This function recalculates iDiffToRef1, diffToRef1, hDistro, and distantSeqs when mutations in
				the loci indicated in toIgnore[] are ignored. In such case, we would assume that the sequences
				are missing all those nucleotides; thus mutations in those positions do not matter. Two
				genotypes will remain connected if they have only one mutation anywhere else. 

				After recompunting these quantities we still need to perform two operations: 
					i) Collapse sequences that are equivalent when ignoring toIgnore[] loci. 
					ii) Connect sequences that are connected when ignoring toIgnore[] loci. 
				These tasks are addressed somewhere else. 

				This function has been tested and debugged to the best of my knowledge! 
				
				Inputs: 
					>> toIgnore: List of nucleotides that we wish to ignore. 

			*/

			int iDummy, iDummy_; 	// These indexes will go over iDiffToRef1[] and toIgnore[] respectively. 
			vector<int> iToErase; 	// This will store the indexes that must be erased in each case. 
			bool fGo=false; 
			seqLength -= toIgnore.size(); 	// Correct sequence length in the renormalized network. 
			// We will recompute distantSeqs from scratch! 
			distantSeqs = new vector<int>[seqLength+1]; 
			for (int iSeq=0; iSeq<seqLength+1; iSeq++){
				distantSeqs[iSeq].clear(); 	// Still absurdly necessary... 
			}

			// Loop over all genotypes -- for each, we look whether mutations w.r.t. ref1 are within toIgnore[]. 
			for (int iSeq=0; iSeq<nGenotypes; iSeq++){
				iDummy=0; 
				iDummy_=0; 
				iToErase.clear(); 
				if (iDiffToRef1[iSeq].size()) fGo = true; // This always is true except for iSeq=refSeq. 
				while (fGo){
					while(iDiffToRef1[iSeq][iDummy]<toIgnore[iDummy_] and iDummy<iDiffToRef1[iSeq].size()) iDummy++; 
					if (iDiffToRef1[iSeq][iDummy]==toIgnore[iDummy_] and iDummy<iDiffToRef1[iSeq].size()){
						iToErase.push_back(iDummy); 
						hDistro[iSeq]--; // This assumes that hDistro[] has been loaded earlier! 
					}
					iDummy_++; 
					if (iDummy_>=toIgnore.size() or iDummy == iDiffToRef1[iSeq].size()) fGo = false; 
				}
				if (hDistro[iSeq]<0) throw std::invalid_argument("ACHTUNG!! Negative number in hDistro[]."); 
				distantSeqs[hDistro[iSeq]].push_back(iSeq); 
				for (int i=(iToErase.size()-1); i>=0; i--){
					iDiffToRef1[iSeq].erase(iDiffToRef1[iSeq].begin()+iToErase[i]); 
					diffToRef1[iSeq].erase(diffToRef1[iSeq].begin()+iToErase[i]); 
				}
			}

			return; 
		}

		void GNet::computeNodeCatalog(){
			/* computeNodeCatalog function: 

				This function computes a catalog of nodes that are grouped with each other. These relationships
				are stored in an object (nodeCatalog[]), that is a vector of vectors. Each vector contains a
				list of nodes that are equivalent (provided that some nucleotides are being ignored). 

				This function assumes that we have already called renormalizeDistros(). It relies on the newly
				computed distantSeqs[], iDiffToRef1[], diffToRef1[], etc, to test whether two sequences became
				the same after ignoring the loci in toIgnore[]. 

				We use distantSeqs[] to compute new IDs. If two sequences have become the same, they are
				necessarily at the same distance w.r.t. reference sequence 0. We only need to compare couples of
				sequences that lie at a same distance then. We loop over all distances to reference 0. The first
				sequence becomes automatically the new nodeID=0. If a sequence has not been matched, it becomes
				a new nodeID in the renormalized network. If a sequence is matched, it becomes associated to the
				corresponding nodeID.

			*/

			// Variables to track whether sequences have already been allocated, and where they belong. 
			nAllGenotypes = nGenotypes; 
			invertedNodeCatalog = new int[nAllGenotypes]; 
			bool* fNew; 
			fNew = new bool[nGenotypes]; 
			for (int iSeq=0; iSeq<nGenotypes; iSeq++){
				fNew[iSeq] = true; 
			}
			nodeCatalog.clear(); 
			int thisSeq1, thisSeq2, iNewID=0; 
			vector<int> dummyVector; 

			// Loop over sequences at a given distance (remind that seqLenth has been recomputed!): 
			for (int iDist=0; iDist<seqLength+1; iDist++){
				// First loop over all sequences at a given distance: 
				for (int iSeq1=0; iSeq1<distantSeqs[iDist].size(); iSeq1++){
					// Has this sequence already been assigned? 
					thisSeq1 = distantSeqs[iDist][iSeq1]; 	// Load this sequence's ID for better handling. 
					if (fNew[thisSeq1]){
						fNew[thisSeq1] = 0; 						// Mark as not new (it will be compared to all others). 
						invertedNodeCatalog[thisSeq1] = iNewID; 	// Mark the renormalized node that it is assigned. 
						nodeCatalog.push_back({thisSeq1});			// Start new vector to store equivalences to thisSeq1. 

						// Second loop over all sequences at a given distance: 
						for (int iSeq2=iSeq1+1; iSeq2<distantSeqs[iDist].size(); iSeq2++){
							thisSeq2 = distantSeqs[iDist][iSeq2]; 	// Load this sequence's ID for better handling. 
							// Has this sequence already been assigned? 
							if (fNew[thisSeq2] and h::areTheSame(iDiffToRef1[thisSeq1], iDiffToRef1[thisSeq2], diffToRef1[thisSeq1], diffToRef1[thisSeq2])){
								fNew[thisSeq2] = false; 
								invertedNodeCatalog[thisSeq2] = iNewID; 
								nodeCatalog[iNewID].push_back(thisSeq2); 
							}
						}
						iNewID++; 	// Increase nodeID for next new node. 
					}
				}
			}

			// Correcting number of genotypes: 
			nGenotypes = nodeCatalog.size(); 

			return; 
		}

		void GNet::computeRepresentativeDistros(){
			/* computeRepresentativeDistros function: 

				This functions computes the relevant distros (hDistro, distantSeqs, iDiffToRef1, diffToRef1)
				assuming that each node is an entry from the nodeCatalog. 

				This also corrects the reference so that now it is a node of the renormalized network. By doing
				this, we are capable of using the functions that we had used before to find connections! 

			*/

			int *dummyHDistro, thisRepSeq; 
			vector<int> *dummyDistantSeqs, *dummyIDiffToRef1; 
			vector<char> *dummyDiffToRef1; 
			dummyHDistro = new int[nGenotypes]; 
			dummyDistantSeqs = new vector<int>[seqLength+1]; 
			for (int i=0; i<seqLength+1; i++){
				dummyDistantSeqs[i].clear(); 
			}
			dummyIDiffToRef1 = new vector<int>[nGenotypes]; 
			dummyDiffToRef1 = new vector<char>[nGenotypes]; 
			for (int iID=0; iID<nGenotypes; iID++){
				thisRepSeq = nodeCatalog[iID][0]; 	// Load a sequence that represents the whole group. 
				dummyHDistro[iID] = hDistro[thisRepSeq]; 
				dummyDistantSeqs[hDistro[thisRepSeq]].push_back(iID); 
				dummyIDiffToRef1[iID].clear(); 
				dummyIDiffToRef1[iID] = iDiffToRef1[thisRepSeq]; 
				dummyDiffToRef1[iID].clear(); 
				dummyDiffToRef1[iID] = diffToRef1[thisRepSeq]; 
			}

			hDistro = dummyHDistro; 
			distantSeqs = dummyDistantSeqs; 
			iDiffToRef1 = dummyIDiffToRef1; 
			diffToRef1 = dummyDiffToRef1; 

			refID1 = invertedNodeCatalog[refID1]; 

			return; 
		}
}

