// Imports:
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <algorithm>
#include <time.h>
using namespace std;

int main(int argc, char *argv[])
{
	srand(time(NULL));
	// Define and initialize variables:
	int n_rw = 1000, tMax_1 = 50000, L = 100, A = 2;
	int iMutate;
	ostringstream fName;
	float epsilon = 0.1;
	fName.str("");
	fName << "./Data/deviation_time.csv";
	ofstream fOut(fName.str());
	fOut << "L"
		 << ","
		 << "mu"
		 << ","
		 << "time"
		 << endl;
	for (int ast = 0; ast < 10; ast++)
	{
		L = (100 + 100 * ast) * 10;
		for (int asta = 0; asta < 10; asta++)
		{
			float mu = 0.1 + (float)asta / 10;
			cout << "Starting L=" << L << " , mu=" << mu << endl;
			// probability of mutation
			int origin = 0; // numer of 1s in origin
			// int tMax = int(tMax_1 * 0.5 * L / (mu * 500));
			int tMax = tMax_1;
			bool display_distros = false;

			int *vocabulary, dummyChar;
			vector<int *> x(n_rw);
			vector<int> reference_sequence(L);
			vocabulary = new int[A];
			vocabulary[0] = 0;
			vocabulary[1] = 1;

			// Initalize the sequence for each random walker
			// x = new char *[n_rw];

			for (int j = 0; j < L; j++)
			{
				if (j < origin)
				{
					reference_sequence[j] = vocabulary[1];
					// cout << reference_sequence[j];
				}
				else
				{
					reference_sequence[j] = vocabulary[0]; // Start sequence (00000...). Modify this if you want a different origin
														   // cout << reference_sequence[j];
				}
			}
			// cout << endl;

			for (int i = 0; i < n_rw; i++)
			{
				x[i] = new int[L];
				for (int j = 0; j < L; j++)
				{
					if (j < origin)
					{
						x[i][j] = vocabulary[1];
					}
					else
					{
						x[i][j] = vocabulary[0]; // Start sequence (00000...). Modify this if you want a different origin
					}
				}
			}
			//   Run the loop:

			float thisD, thisMeanD2, thisD_eucl;
			vector<float> d2, d2_eucl;
			vector<int> desired_t_vec;

			for (int t = 0; t < tMax; t++)
			{

				// cout << t << " of " << tMax << endl;
				d2.clear();
				d2_eucl.clear();
				// cout << "Time: " << t << endl;
				thisMeanD2 = 0;
				vector<vector<int>> aux_x; // This is an auxiliary vector that will contain size-free replicating population

				// For each walker in the existing population
				for (int i = 0; i < n_rw; i++)
				{

					// Create n_copies copies with mutation. Add them to aux_x population
					int n_copies = 1; // This will be changed to 2 in some cases in the future.
					for (int k = 0; k < n_copies; k++)
					{
						// Select the i-th random walker
						vector<int> rw(L);
						for (int j = 0; j < L; j++)
						{
							rw[j] = x[i][j];
						}

						// Create single-point mutation on the sequence with probability mu
						float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
						if (r <= mu)
						{
							iMutate = rand() % L;

							rw[iMutate] = 1 - rw[iMutate];
						}

						// Save this sequence in auxiliar population
						aux_x.push_back(rw);
					}
				}

				// Re-sample the population to stay with constant size n_rw
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

				// Compute distances
				// 1- Hamming distance (d_H)
				// 2- Euclidean distance   (d_E)
				// 3- Ortogonal-vectors distance (d_bi)

				for (int i = 0; i < n_rw; i++)
				{
					thisD = 0;
					thisD_eucl = 0;
					for (int j = 0; j < L; j++)
					{
						// cout << x[i][j] << endl;
						thisD += fabs(x[i][j] - reference_sequence[j]);
						thisD_eucl += pow(fabs(x[i][j] - reference_sequence[j]), 2);
					}
					// thisMeanD2 += thisD; // pow(thisD, 2);
					d2.push_back(thisD); // Save the distance-squared value in d2
					d2_eucl.push_back(sqrt(thisD_eucl));
				}

				// Compute the mean distance

				float mean_d = 0, mean_d_eucl = 0;
				for (int d_index = 0; d_index < n_rw; d_index++)
				{
					float d_eucl = d2[d_index];
					mean_d += d_eucl;
					mean_d_eucl += float(d2_eucl[d_index]);

					// cout << d_eucl << endl;
				}
				mean_d = mean_d / n_rw;
				mean_d_eucl = mean_d_eucl / n_rw;

				if (mean_d <= int(mu * (1 - epsilon) * float(t + 1)))
				{
					fOut << L << "," << mu << "," << t + 1 << endl;
					cout << "Breaking at t=" << t + 1 << endl;
					break;
				}
			}

			delete[] vocabulary;
			for (int i = 0; i < n_rw; i++)
			{
				delete[] x[i];
			}
		}
	}

	return 0;
}
