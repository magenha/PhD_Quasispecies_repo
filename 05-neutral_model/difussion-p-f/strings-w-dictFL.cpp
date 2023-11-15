// Imports:
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <time.h>
#include <tuple>

using namespace std;

int Humming_distance(vector<char> s1, vector<char> s_ref)
{
	int thisD = 0;
	int L = s1.size();
	for (int j = 0; j < L; j++)
	{
		// //cout << x[i][j] << endl;
		if (s1[j] != s_ref[j])
			thisD += 1;
	}

	return thisD;
}

int main(int argc, char *argv[])
{

	srand(time(NULL));
	// Define and initialize variables:
	int n_rw = 1000, L = 10000, A = 4;
	int tMax = 20; // other option: int(0.005 * float(L))
	int iMutate;   // This is a dummy variable to indicate the locus of the mutation

	float p = 1.0; // Holeyness of the FL
	vector<int> f_list;
	f_list.push_back(1);
	f_list.push_back(2);
	f_list.push_back(4);
	f_list.push_back(8);
	for (int zz = 0; zz < f_list.size(); zz++)
	{
		cout << "Fitness = " << f_list[zz] << endl;
		int f = f_list[zz];						   // This is the value of the fitness when no nule value is used
		vector<tuple<vector<char>, int>> f_values; // This is the dict used to store each sequence's fitness
		char *vocabulary, dummyChar;
		vector<char *> x(n_rw); // This vector is the population
		vocabulary = new char[A];
		vocabulary[0] = 'A';
		vocabulary[1] = 'C';
		vocabulary[2] = 'G';
		vocabulary[3] = 'T';

		// Define the Output file name
		ostringstream fName;
		fName.str("");
		fName
			<< "./Data/DictFL/DictFL_p_" << p << "_f_" << f << "_.csv";
		ofstream fOut(fName.str());

		// Initalize the sequence for each random walker
		for (int i = 0; i < n_rw; i++)
		{
			// cout << "Initiating sequence " << i << ": ";
			x[i] = new char[L];
			for (int j = 0; j < L; j++)
			{
				x[i][j] = vocabulary[0];
				// cout << vocabulary[0];
			}
			// cout << endl;
		}
		// cout << "Adding f to wild type sequence " << endl;
		// cout << "To check: ";
		//  Set the wt with f descendants
		vector<char> sequence(L);
		for (int m = 0; m < L; m++)
		{
			sequence[m] = x[0][m];
			// cout << sequence[m];
		}
		tuple<vector<char>, int> new_element(sequence, f);
		f_values.push_back(new_element);
		// cout << " Descendants=" << get<1>(f_values[0]) << endl;

		// Run the loop:
		float thisD,
			thisMeanD2; // ThisD computes the distance of the sequence. ThisMeanD2 cummulates the sum for al random walkers.
		vector<float> d2;
		vector<vector<char>> aux_x; // This is an auxiliary vector that will contain size-free replicating population

		// In each step of time
		for (int t = 0; t < tMax; t++)
		{
			cout << "time step=" << t + 1 << endl;
			// cout << "Size of fitness dict " << f_values.size() << endl;
			d2.clear();
			aux_x.clear();
			thisMeanD2 = 0;

			// cout	<< "size of auxiliar distances values vector " << d2.size() << " ";
			// cout << "size of the auxiliar population vector " << aux_x.size() << endl;

			// For each walker in the existing population
			for (int i = 0; i < n_rw; i++)
			{
				// Detect number of descendants for that sequence
				// Check if the sequence is in f_values
				// cout << "walker " << i << " " << endl;
				int n_desc;		   // variable to know the assigned number of descandants
				bool assigned = 0; // variable to know if the sequence exists or not in our dict
				// cout << "Searching if it is in dict" << endl;
				//  Search if the sequence is in dict
				for (int j = 0; j < f_values.size(); j++)
				{
					// auto seq = f_values[j];

					// cout << "Dict sequence number " << j << endl;

					bool all_equal = true; // Flag to check if all characters are equal, if true then sequence is in dict, false means that is not in dict
					// iterate over all chars of the sequence
					for (int k = 0; k < L; k++)
					{
						if (x[i][k] != get<0>(f_values[j])[k])
						{
							all_equal = false;
							break; // Exit the loop if any character is not equal
						}
					}

					if (all_equal) // sequence in dict
					{
						n_desc = get<1>(f_values[j]);
						assigned = 1;
						break;
					}
				}

				// If the sequence is not in dict, then add it
				if (assigned == 0)
				{
					// cout << "Sequence not in dict";
					//  Compute the assigned replicability
					float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
					if (r <= p)
					{
						// cout << "r value = " << r << ", then positive fitness" << endl;
						n_desc = f;
						vector<char> sequence(L);
						for (int m = 0; m < L; m++)
						{
							sequence[m] = x[i][m];
						}
						// Pass the x[i] value to a vector<char *>
						tuple<vector<char>, int> new_element(sequence, n_desc);
						f_values.push_back(new_element);
					}
					else
					{
						// cout << "r value = " << r << ", then zero fitness" << endl;
						n_desc = 0;
						vector<char> sequence(L);
						for (int m = 0; m < L; m++)
						{
							sequence[m] = x[i][m];
						}
						// Pass the x[i] value to a vector<char *>
						tuple<vector<char>, int> new_element(sequence, n_desc);
						f_values.push_back(new_element);
					}
				}

				else
				{
					// cout << "Sequence already in dict, f=" << n_desc << endl;
				}

				// Create f copies with mutation. Add them to aux_x population
				if (n_desc > 0)
				{

					for (int k = 0; k < n_desc; k++)
					{
						// Copy i-th random walker
						vector<char> rw(L);
						for (int j = 0; j < L; j++)
						{
							rw[j] = x[i][j];
						}
						// Create single-point mutation on that
						iMutate = rand() % L;
						dummyChar = vocabulary[rand() % A];
						while (dummyChar == x[i][iMutate]) // Make sure you do not select the same nucleotide
						{
							dummyChar = vocabulary[rand() % A];
						}
						rw[iMutate] = dummyChar;

						// Save this sequence in auxiliar population
						aux_x.push_back(rw);
					}
				}
				else
				{
					// cout << "No descendants, then no mutations to create, just appending" << endl;
					vector<char> rw(L);
					for (int j = 0; j < L; j++)
					{
						rw[j] = x[i][j];
					}
					aux_x.push_back(rw);
				}
			}

			// This might be a sepparated function
			// Now select just n_rw random sequences within aux_x
			int sizeofAuxX = aux_x.size();
			// cout << "Size of aux_X " << sizeofAuxX << endl;
			for (int i = 0; i < n_rw; i++)
			{
				int irw = rand() % sizeofAuxX;
				for (int j = 0; j < L; j++)
				{
					x[i][j] = aux_x[irw][j];
				}
			}

			// This might be a sepparate function & section
			// devoted to analyze the population in the given moment
			// of the process

			// compute average number of mutations
			for (int i = 0; i < n_rw; i++)
			{
				thisD = 0;
				for (int j = 0; j < L; j++)
				{
					// //cout << x[i][j] << endl;
					if (x[i][j] != vocabulary[0])
						thisD += 1;
				}
				thisMeanD2 += thisD;
				d2.push_back(thisD); // This needs to be stored for the variance
			}

			// This is temporal, a check for the d2 vector
			// ostringstream fName2;
			// fName2.str("");
			// fName2
			//	<< "./Data/test_d2_" << t << "_.csv";
			// ofstream fOut2(fName2.str());
			// for (int o = 0; o < d2.size(); o++)
			//{
			//	fOut2 << d2[o] << endl;
			//}

			thisMeanD2 = thisMeanD2 / n_rw;
			// //cout << thisMeanD2 << endl;

			// Compute the variance
			float variance_d = 0;
			for (int d_index = 0; d_index < d2.size(); d_index++)
			{
				variance_d += pow(d2[d_index] - thisMeanD2, 2);
			}
			variance_d = variance_d / d2.size();

			//
			//  Maybe save the values of <m> and Var(m) ?
			fOut << thisMeanD2 << "," << variance_d << endl;
			//

			// cout << "\n\n\n\n"<< endl;
		}

		// Check that fraction of null vs. no null values in f is p.
		float f_sum = 0;
		for (int i = 0; i < f_values.size(); i++)
		{
			f_sum += get<1>(f_values[i]);
		}
		cout << "Fitness Dict has mean value " << f_sum / (float(f_values.size())) << " . And p=" << p << endl;

		// Free some of the memory usage
		// Check this properly
		delete[] vocabulary;
		for (int i = 0; i < n_rw; i++)
		{
			delete[] x[i];
		}
	}

	return 0;
}
