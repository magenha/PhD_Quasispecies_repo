

// Imports:
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <sstream>
using namespace std;

int main(int argc, char *argv[])
{

	// Define and initialize variables:
	int n_rw = 1000, tMax = 10000, L = 50000;
	int iMutate;
	char *vocabulary, **x, dummyChar;
	vocabulary = new char[4];
	vocabulary[0] = 'A';
	vocabulary[1] = 'C';
	vocabulary[2] = 'G';
	vocabulary[3] = 'T';
	x = new char *[n_rw];
	for (int i = 0; i < n_rw; i++)
	{
		x[i] = new char[L];
		for (int j = 0; j < L; j++)
		{
			x[i][j] = 'A';
		}
	}
	ostringstream fName;
	fName.str("");
	fName << "./dataOut.csv";
	ofstream fOut(fName.str());

	// Run the loop:
	float thisD, thisMeanD2;
	vector<float> d2;
	d2.clear();
	for (int t = 0; t < tMax; t++)
	{
		thisMeanD2 = 0;
		for (int i = 0; i < n_rw; i++)
		{
			iMutate = rand() % L;
			dummyChar = vocabulary[rand() % 4];
			while (dummyChar == x[i][iMutate])
			{
				dummyChar = vocabulary[rand() % 4];
			}
			x[i][iMutate] = dummyChar;
			thisD = 0;
			for (int j = 0; j < L; j++)
			{
				if (x[i][j] != 'A')
					thisD += 1;
			}
			thisMeanD2 += thisD; // pow(thisD, 2);
		}
		cout << thisMeanD2 / n_rw << endl;
		fOut << thisMeanD2 / n_rw << endl;
		d2.push_back(thisMeanD2 / n_rw);
	}

	fOut.close();

	return 0;
}
