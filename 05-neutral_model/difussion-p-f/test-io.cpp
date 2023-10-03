

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
	vocabulary[0] = 'C';
	vocabulary[0] = 'G';
	vocabulary[0] = 'T';
	x = new char *[n_rw];
	for (int i = 0; i < n_rw; i++)
	{
		x[i] = new char[L];
		for (int j = 0; j < L; j++)
		{
			x[i][j] = 'A';
		}
	}

	//// For I-O:

	// Las variables de tipo ostringstream son una generalización de los strings.
	// Permiten hacer muchas operaciones de manera sencilla.
	ostringstream fName;
	fName.str(""); // Esto limpia el contenido de la variable, ya que podría estar almancenando algo en la posición de memoria.
	fName << "./dataOut.csv";

	// Las variables de tipo ofstream son las que se utilizan para hacer output de datos.
	// Para hacer input, se usan variables de tipo ifstream, pero de eso no pongo nada en este código.
	ofstream fOut(fName.c_str()); // Esto crea la variable fOut (de tipo ofstream) y abre el archivo con el nombre guardado en fName.

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
		fOut << thisMeanD2 / n_rw << endl; // Ahora puedes tratar la variable fOut como si fuera cOut.
	}

	fOut.close(); // Cerramos el archivo para que el compilador no se estrese.

	return 0;
}
