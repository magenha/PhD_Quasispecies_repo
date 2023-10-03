This folder contains code two reconstruct genotype-genotype networks, with two examples.

In order for the code to work properly, the following steps are needed:

	1) Let ~/ be the home folder where all code resides.
	2) Within ~/Data/, we need to have one .fasta file for each sample of genotypes.
		>> In the examples provided, there are 5 files (one from each major variant of the SARS-CoV-2 virues).
		>> Each of these files contains haplotypes of sequences in a specific viral region.
	3) Within ~/Data/, we need to create a file called fileNames.csv that contains the names (without extensions) of each .fasta file.
		>> Make sure that there isn't an empty line at the end of the file.
	4) Now we can run "script_findIDs.py":
		>> This script finds unique sequences across all files, and produces a bare list of sequences.
		>> These sequences are listed in the file "nodesID.csv".
		>> They appear in the order in which they were found.
		>> Each node's ID is just the order (starting at 0) in which they appear in this file.
	5) Now we can run the script "script_extractAbundanceID.py":
		>> This script creates a folder ~/Data/Abundances/ where, for each variant, it stores:
			> nodeID, abundance (as absolute number), abundance (as frequency), abundance (as component of a unitary vector).
		>> The additional forms of abundance are useful, e.g., to compute Kullback-Leibler divergences or cosine distances.
	6) Within  ~/Data/, now we need to create a file called "metadata.csv":
		>> There is an example of this file in each folder.
		>> It must contain three lines:
			> First, containing just the number of genotypes in this region (i.e. number of lines in "nodesID.csv").
			> Second, containing the sequence that will be used as a reference to triangulate the nodes.
				- ID 0 works fine. For larger network, it is advisable to find a central node (usually one with large abundance).
			> Third, containing the length of the sequence (number of characters in each line of "nodesID.csv").
	7) Now we can compile and run the script "cScript_generateNetwork_.cpp":
		>> In line 60 of this script there is a variable netName:
			> This is not relevant, but we might consider changing it to something informative.
			> In the examples provided, each region is correctly identified, so that the script can be compiled and run rightaway.
		>> This script generates several important pieces of information:
			> A folder ~/Data/GN/ to contain all the info that the script outputs.
			> A file "diffsToRef_0.csv" with the mutations of each sequence w.r.t. the first reference (in the examples, ID 0).
				- ACHTUNG!! This file is written in a specific format for the scripts. It is not readily understandable.
			> A file "distantDistro.csv" with the distributions of distances to the first reference, with the following structure:
				- First, sequence length (which marks the longest possible distance to the sequence).
				- Second, number of sequences at distance 0 (usually just the reference).
				- Third, list of nodes at distance 0, listed each in one line (in case of distance 0, just the ID of the reference sequence).
				- Fourth, number of sequences at distance 1 (probably more than one).
				- Fifth, list of nodes at distance 1, listed each in one line.
				- etc etc etc.
			> A file "hDistro_ref1.csv" that contains the same information as "distantDistro.csv" written differently:
				- Now, each line corresponds to a sequence and it contains the distance of that sequence to the first reference.
			> A similar file "hDistro_ref2.csv" that contains this information for the second reference:
				- The second reference is the sequence furthest away from the first reference, and it is computed internally by the script.
			> A file "edges.csv" containing a list of edges in the network, with nodes identified with their ID as usual:
				- Note that isolated nodes will not appear here, since they do not participate of any edge.
			> A file "grandNetwork.csv" with the network stored in a format that is specific for the network reconstruction software:
				- This file contains the same information as "edges.csv" and, additionally, it loads nodes that are isolated.
			> A file "whereMutations.csv" that stores where, along the sequence happens each mutation of connected nodes:
				- This information is stored in a fairly complicated way.
				- The only way to access this information at the moment is to load "grandNetwork.csv" as a C++ object as programed in the network reconstruction software. Only then does "whereMutations.csv" make sense.
				- On the happy side, there are functions implemented to load this information in the gNet.h C++ library.
	8) Now we can run the script "script_generateSubNetworks.py", that extracts the subnetwork corresponding to each variant:
		>> An important variable is nLoad:
			> This controls whether we want to built the complete subnetwork for each variants, or just subnetworks with the top nLoad most abundance sequences of each variant.
			> If nLoad<0:
				- It extracts complete subnetworks.
			> If nLoad>0:
				- It extracts subnetworks with just the top nLoad most abundant sequences of each variante.
				- It then proceeds to generate an aggregate network of these subnetworks (note this is smaller than the grandNetwork).
		>> This generates a folder ~/Data/GN/SubNetworks/ that contains one file per variant with the edges of the corresponding subnetwork.
