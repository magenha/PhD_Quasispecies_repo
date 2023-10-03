/*

	helper.h: 

		This file contains the library helper to aid with the calculations of genotype networks. 

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
using namespace std; 
// std::default_random_engine kMapGenerator; 


namespace h{

	/*

		helper namespace. 

	*/ 

	int getIntFromString(string toRead){
		/* getIntFromString function: 

			This function returns the integer stored in the provided string. 

		Inputs: 
			>> toRead: string to be converted into an integer. 

		Returns: 
			>> intValue: integer value stored in toRead. 

		*/ 

		int intValue;
		stringstream sLine;

		sLine.clear();
		sLine.str(toRead);
		sLine >> intValue;

		return intValue;
	}

	float getFloatFromString(string toRead){
		/* getIntFromString function: 

			This function returns the integer stored in the provided string. The difference with the
			previous function is that it returns a float variable. 

		Inputs: 
			>> toRead: string to be converted into an integer. 

		Returns: 
			>> intValue: integer value stored in toRead. 

		*/ 

		float floatValue;
		stringstream sLine;

		sLine.clear();
		sLine.str(toRead);
		sLine >> floatValue;

		return floatValue;
	}

	int hammingDistance(int seqLength, string seq1, string seq2){
		/* hammingDistance function: 

			This function computes the Hamming distance between two sequences. 

			Inputs: 
				>> seqLength: length of the sequences provided. 
				>> seq1, seq2: Of which we will measure the distance. 

			Returns: 
				<< hDistance: Hamming distance between sequences seq1 and seq2. 

		*/ 

		int hDistance=0; 
		for (int i=0; i<seqLength; i++){
			if (seq1[i] != seq2[i]) hDistance++; 
		}

		return hDistance; 
	}

	int hammingDistanceFromDifferences(vector<int> iDiff1, vector<int> iDiff2, vector<char> diff1, vector<char> diff2){
		/* hammingDistanceFromDifferences function: 

			This function computes the Hamming distance between two strings given the differences of each string to a reference
			string. This function is thought out to speed up the computation of Hamming distances to a second reference string
			once distances and differences of all strings to a first reference have already been computed. 

			Inputs: 
				>> iDiff1, iDiff1: Vector of integers indicating the positions at which each sequence differs from the first
				reference. 
				>> diff1, diff2: Vector with the nucleotides at each position where each sequence differs from the first reference. 

			Returns: 
				<< hDistance: Hamming distance between both sequences. 

		*/ 

		int hDistance=0, k1=0, k2=0; 
		float fGo=true; 

		// Loop over all entries in iDiff1 and iDiff2: 
		while(fGo){
			while (iDiff1[k1] != iDiff2[k2] and ((k1 < iDiff1.size()) and (k2 < iDiff2.size()))){
				while ((iDiff1[k1] < iDiff2[k2]) and (k1 < iDiff1.size())){
					k1++; 
					hDistance++; 
				}
				while ((iDiff1[k1] > iDiff2[k2]) and (k2 < iDiff2.size())){
					k2++; 
					hDistance++; 
				}
			}
			if ((k1 < iDiff1.size()) and (k2 < iDiff2.size())){
				if (diff1[k1] != diff2[k2]) hDistance++; 
				k1++; 
				k2++; 
			}
			if ((k1 == iDiff1.size()) or (k2 == iDiff2.size())) fGo=false; 
		}
		if (k1<iDiff1.size()) hDistance += iDiff1.size()-k1; 
		if (k2<iDiff2.size()) hDistance += iDiff2.size()-k2; 

		return hDistance; 
	}

	vector<int> hammingDistanceAndDifferences(int seqLength, string seq1, string seq2){
		/* hammingDistanceAndDifferences functino: 

			This function computes the Hamming distance between two sequences and the positions along the sequence where the
			differences happen. All this information is passed along in the returned vector of integers: entries 0:(hDistance -
			1) contain the positions at which seq1 and seq2 differ. Position hDistance constains the Hamming distance. This
			information has to be unpacked outside of the function. 

			Inputs: 
				>> seqLength: length of the sequences provided. 
				>> seq1, seq2: Of which we will measure the distance. 

			Returns: 
				<< summaryVector: A vector containing the positions at which both sequences differ and their Hamming distance. 

		*/ 

		vector<int> summaryVector; 
		summaryVector.clear(); 
		int hDistance=0; 
		for (int i=0; i<seqLength; i++){
			if (seq1[i] != seq2[i]){
				hDistance++; 
				summaryVector.push_back(i); 
			}
		}
		summaryVector.push_back(hDistance); 

		return summaryVector; 
	}

	bool areConnected(int seqLength, string seq1, string seq2, int hThreshold=1){
		/* areConnected function: 

			This function runs over sequences seq1 and seq2 and returns True if they are not further apart than the provided
			threshold. 

			Inputs: 
				>> seqLength: length of the provided sequences. 
				>> seq1, seq2: Sequences to be compared. 
				>> hThreshold=1: Maximum distance that considers seq1 and seq2 connected. 

			Returns: 
				<< true, if h(seq1, seq2) <= hThreshold; false otherwise. 

		*/ 

		int h = 0; 
		for (int i=0; i<seqLength; i++){
			if (seq1[i] != seq2[i]) h++; 
			if (h > hThreshold) return false; 
		}

		return true; 
	}

	bool areConnected(vector<int> iDiff1, vector<int> iDiff2, vector<char> diff1, vector<char> diff2){
		/* areConnected function: 

			This function determines whether two sequences are connected by looking only at the places at which they differ from
			a reference sequence. 

			Inputs: 
				>> iDiff1, iDiff2: Vector of integers with the indexes of the locations in which each sequence differs from the
				reference sequence. 
				>> diff1, diff2: Vector of characters with the nucleotides at the locations in which each sequence differs from the
				reference sequence. 

			Returns: 
				<< true, if both sequences differ in no more than one nucleotide; false otherwise. 

		*/ 

		// Auxiliary variable to track state: 
		bool auxFlag=false; 

		// If both sequences differ in the same amount of positions: 
		if (iDiff1.size() == iDiff2.size()){
			// We need to check that the positions in which they differ are exactly the same, and that only one nucleotide is
			// different. 
			for (int i=0; i<iDiff1.size(); i++){
				// If differences are not in all same nucleotides, there are at least 2 differences, sequences are not connected. 
				if (iDiff1[i] != iDiff2[i]) return false; 
				// If differences are in same nucleotide, but nucleotide differs: 
				if (diff1[i] != diff2[i]){
					// This can only happen once! If it happened before, sequences are not connected. 
					if (auxFlag) return false; 
					auxFlag = true; 
				}
			}
		}
		// If both sequences have different lengths: 
		else{
			// Dummy variables to correct for size difference of iDiff and diff vectors: 
			int minDiff, k1=0, k2=0; 
			minDiff = min(iDiff1.size(), iDiff2.size()); 
			// Loop over differing nucleotides: 
			for (int i=0; i<minDiff; i++){
				// Because of different sizes, there will be at most one case in which connected sequences present a different site: 
				if (iDiff1[i+k1] != iDiff2[i+k2]){
					// If this has happened before, sequences are not connected! 
					if (auxFlag) return false; 
					auxFlag = true; 
					// If seq1 has more mutations than seq 2, we resume comparison a position ahead -- and vice versa: 
					if (iDiff1.size()>iDiff2.size()) k1=1; 
					else k2=1; 
					// If mutation sites still differ after correction, sequences are not connected: 
					if (iDiff1[i+k1]!=iDiff2[i+k2]) return false; 
				}
				// The only difference must be on the extra site. All other mutations must be the same. Otherwise, not connected: 
				if (diff1[i+k1] != diff2[i+k2]) return false; 
			}
		}

		return true; 
	}

	int areConnected_(vector<int> iDiff1, vector<int> iDiff2, vector<char> diff1, vector<char> diff2){
		/* areConnected function: 

			This function determines whether two sequences are connected by looking only at the places at which they differ from
			a reference sequence. 

			What follows is just magic: 
				Besides returning whether two sequences are connected or not, I now need to return the position at which the
				mutation happens. But c++ is notably bad at allowing to return more than one value. Therefore, if two sequences are
				connected because they have a mutation at position i, this function will return numerical value i+1. This will be
				interpreted as true -- even for connections due to a mutation at position 0. Thus we can recover all the information
				needed. 

			Inputs: 
				>> iDiff1, iDiff2: Vector of integers with the indexes of the locations in which each sequence differs from the
				reference sequence. 
				>> diff1, diff2: Vector of characters with the nucleotides at the locations in which each sequence differs from the
				reference sequence. 

			Returns: 
				<< (i+1), where i is the position of the mutation if both sequences differ in no more than one nucleotide. 
				<< 0, otherwise. 

		*/ 

		// Auxiliary variable to track state: 
		bool auxFlag=false; 
		int iMu = -1; 

		// If both sequences differ in the same amount of positions: 
		if (iDiff1.size() == iDiff2.size()){
			// We need to check that the positions in which they differ are exactly the same, and that only one nucleotide is
			// different. 
			for (int i=0; i<iDiff1.size(); i++){
				// If differences are not in all same nucleotides, there are at least 2 differences, sequences are not connected. 
				if (iDiff1[i] != iDiff2[i]) return 0; 
				// If differences are in same nucleotide, but nucleotide differs: 
				if (diff1[i] != diff2[i]){
					// This can only happen once! If it happened before, sequences are not connected. 
					if (auxFlag) return 0; 
					auxFlag = true; 
					iMu = iDiff1[i]; 
				}
			}
		}
		// If both sequences have different lengths: 
		else{
			// Dummy variables to correct for size difference of iDiff and diff vectors: 
			int minDiff, k1=0, k2=0; 
			minDiff = min(iDiff1.size(), iDiff2.size()); 
			// Loop over differing nucleotides: 
			for (int i=0; i<minDiff; i++){
				// Because of different sizes, there will be at most one case in which connected sequences present a different site: 
				if (iDiff1[i+k1] != iDiff2[i+k2]){
					// If this has happened before, sequences are not connected! 
					if (auxFlag) return 0; 
					auxFlag = true; 
					// If seq1 has more mutations than seq 2, we resume comparison a position ahead -- and vice versa: 
					if (iDiff1.size()>iDiff2.size()){
						iMu = iDiff1[i]; 
						k1=1; 
					}
					else{
						iMu = iDiff2[i]; 
						k2=1; 
					}
					// If mutation sites still differ after correction, sequences are not connected: 
					if (iDiff1[i+k1]!=iDiff2[i+k2]) return 0; 
				}
				// The only difference must be on the extra site. All other mutations must be the same. Otherwise, not connected: 
				if (diff1[i+k1] != diff2[i+k2]) return 0; 
			}
			// Nice bug here!! 
			// 	If the diferring nucleotide is the last one in the sequence that has more mutations, then the above loop is
			// 	finished and it is clear that both seqs. are connected. However, because we never looped until the end of the
			// 	iDiff for the seq. with more mutations, we never read the position of that mutation. 
			// Fix:
			if (iMu<0){
				if (iDiff1.size()>iDiff2.size()) iMu = iDiff1[iDiff1.size()-1]; 
				else iMu = iDiff2[iDiff2.size()-1]; 
			}
		}

		return iMu+1; 
	}

	bool areTheSame(vector<int> iDiff1, vector<int> iDiff2, vector<char> diff1, vector<char> diff2){
		/* areTheSame function: 

			This function computes whether two sequences are the same based on the differences of each of
			them w.r.t. a suitable reference. 

			This function is much simpler than areConnected() because we just need to make sure that all
			mutations happen at exactly the same loci, and that all mutations changed into the same
			nucleotide. 

			Inputs: 
				>> iDiff1, iDiff2: Vector of integers with the indexes of the locations in which each sequence differs from the
				reference sequence. 
				>> diff1, diff2: Vector of characters with the nucleotides at the locations in which each sequence differs from the
				reference sequence. 

			Returns: 
				<< false: If there is any difference at all between iDiff1 and iDiff2, or between diff1 and
				diff2. 
				<< true: If both iDiff1 and iDiff2 are the same and diff1 and diff2 are the same. 
				
		*/ 

		for (int i=0; i<iDiff1.size(); i++){
			if (iDiff1[i] != iDiff2[i]) return false; 
			if (diff1[i] != diff2[i]) return false; 
		}

		return true; 
	}

	string loadGenotype(int gID, string fInName){
		/* loadGenotype function: 

			This function loads genotype number gID from the provided file. Therefore, it makes use of grep and command line. 

			Inputs: 
				>> gID: ID of the genotype that we want to retrieve. 
				>> fInName: File containing the sorted genotypes. 

			Returns: 
				<< seq: Sequence of the retrieved genotype. 

		*/ 

		// Extracting sequence to temporal file: 
		ostringstream strCall; 
		strCall.str(""); 
		strCall << "sed '" << gID+1 << "!d' "  << fInName.c_str() << " > temp.txt"; 
		system(strCall.str().c_str()); 

		// Reading sequence from temporal file: 
		strCall.str(""); 
		strCall << "temp.txt"; 
		ifstream fIn(strCall.str().c_str()); 
		string seq; 
		getline(fIn, seq); 

		return seq; 
	}

	vector<int> rProtocol_nConsecutive(int nMiss, int seqLength, string fOutBaseName="./dummyOut"){
		/* rProtocol_nConsecutive function: 

			This function generates a renormalization protocol in which nMiss consecutive nucleotides are
			missing. 

			The resulting protocol is stored in two files: fOutBaseName+"_keep.csv" and
			fOutBaseName+"_erase.csv". This is for the record. 

			Inputs: 
				>> nMiss: number of nucleotides that will be omitted. 
				>> seqLength: length of the sequence that will be renormalized. 
				>> fOutBaseName: base name to store the generated protocol. 

			Returns: 
				<< iErase: vector with the positions of the nucleotides that will be ignored. 

		*/

		// Files where the protocol will be stored: 
		ostringstream fOutName; 
		fOutName.str(""); 
		fOutName << fOutBaseName << "_keep.csv"; 
		ofstream fOutKeep(fOutName.str().c_str()); 
		fOutName.str(""); 
		fOutName << fOutBaseName << "_erase.csv"; 
		ofstream fOutErase(fOutName.str().c_str()); 

		vector<int> iErase; 
		iErase.clear(); 
		int iStart=rand()%(seqLength-nMiss+1); 
		
		int i; 
		for (i=0; i<iStart; i++){
			fOutKeep << i << endl; 
		}
		for (i; i<iStart+nMiss; i++){
			fOutErase << i << endl; 
			iErase.push_back(i); 
		}
		for (i; i<seqLength; i++){
			fOutKeep << i << endl; 
		}

		return iErase; 
	}

}


