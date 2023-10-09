// Imports:
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <time.h>
using namespace std;

int main(int argc, char *argv[])
{
	srand(time(NULL));
	// Define and initialize variables:
	int n_rw = 1000, tMax = 10000, L = 1000, A = 4;
	int iMutate;
	int f = 1; // replication
	char *vocabulary, dummyChar;
	vector<char *> x(n_rw);
	vocabulary = new char[A];
	vocabulary[0] = 'A';
	vocabulary[1] = 'C';
	vocabulary[2] = 'G';
	vocabulary[3] = 'T';

	ostringstream fName;
	fName.str("");
	string fstr = to_string(f);
	fName
		<< "./dataOut_f_" << fstr << ".csv";
	ofstream fOut(fName.str());

	// Initalize the sequence for each random walker
	// x = new char *[n_rw];
	for (int i = 0; i < n_rw; i++)
	{
		x[i] = new char[L];
		for (int j = 0; j < L; j++)
		{
			x[i][j] = vocabulary[0];
		}
	}
	//   Run the loop:
	float thisD, thisMeanD2;
	vector<float> d2;

	for (int t = 0; t < tMax; t++)
	{
		cout << t << endl;
		d2.clear();
		// cout << "Time: " << t << endl;
		thisMeanD2 = 0;
		vector<vector<char>> aux_x; // This is an auxiliary vector that will contain size-free replicating population

		// For each walker in the existing population
		for (int i = 0; i < n_rw; i++)
		{
			// Create f copies with mutation. Add them to aux_x population
			for (int k = 0; k < f; k++)
			{
				// Select the i-th random walker
				vector<char> rw(L);
				for (int j = 0; j < L; j++)
				{
					rw[j] = x[i][j];
				}
				// Create single-point mutation
				iMutate = rand() % L;
				dummyChar = vocabulary[rand() % A];
				while (dummyChar == x[i][iMutate])
				{
					dummyChar = vocabulary[rand() % A];
				}
				rw[iMutate] = dummyChar;

				// Save this sequence in auxiliar population
				aux_x.push_back(rw);
			}
		}
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

		// compute distances
		for (int i = 0; i < n_rw; i++)
		{
			thisD = 0;
			for (int j = 0; j < L; j++)
			{
				// cout << x[i][j] << endl;
				if (x[i][j] != vocabulary[0])
					thisD += 1;
			}
			thisMeanD2 += thisD; // pow(thisD, 2);
			d2.push_back(thisD); // Save the distance-squared value in d2
		}

		// Compute the mean distance
		float mean_d = 0;
		for (int d_index = 0; d_index < n_rw; d_index++)
		{
			float d_eucl = sqrt(d2[d_index]);
			mean_d += d_eucl;
		}
		mean_d = mean_d / n_rw;

		// Compute the variance
		float variance_d = 0;
		for (int d_index = 0; d_index < d2.size(); d_index++)
		{
			float d_eucl = sqrt(d2[d_index]);
			variance_d += pow(d_eucl - mean_d, 2);
		}
		variance_d = variance_d / d2.size();

		// cout << mean_d << "," << variance_d << "," << thisMeanD2 / n_rw << endl;
		fOut << mean_d << "," << variance_d << "," << thisMeanD2 / n_rw << endl;
	}

	delete[] vocabulary;
	for (int i = 0; i < n_rw; i++)
	{
		delete[] x[i];
	}

	return 0;
}
