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

int main(int argc, char *argv[])
{
    srand(time(NULL));
    vector<char> seq(2);
    vector<tuple<vector<char>, int>> tup_vec(10);
    seq[0] = 'A';
    seq[1] = 'A';
    char *vocabulary, dummyChar;
    vocabulary = new char[4];
    vocabulary[0] = 'A';
    vocabulary[1] = 'C';
    vocabulary[2] = 'G';
    vocabulary[3] = 'T';
    for (int j = 0; j < 10; j++)
    {
        seq[random() % 2] = vocabulary[random() % 4];
        tup_vec[j] = make_tuple(seq, j);
        auto tup = tup_vec[j];
        auto a = get<0>(tup);
        for (int i = 0; i < 2; i++)
        {
            cout << a[i];
        }
        cout << ", " << get<1>(tup) << endl;
    }

    vector<char> desired(2);
    desired[0] = 'A';
    desired[1] = 'A';

    return 0;
}

int is_in_dict(vector<char> desired, vector<tuple<vector<char>, int>> tup_vec)
{
    /*
        This function evaluates if a given sequence is in a given
        fitness dict or not.
        Returns
            1 ->  if is in dict
            0 ->  if is not in dict
    */
    int dummyFlag = 0;
    int s = tup_vec.size();
    int s_des = desired.size();
    for (int j = 0; j < s; j++)
    {
        dummyFlag = 0;
        for (int k = 0; k < s_des; k++)
        {
            auto seq = get<0>(tup_vec[j]);
            if (desired[k] == seq[k])
            {
                dummyFlag += 1;
            }
        }
        if (dummyFlag == s_des)
        {
            return 1;
        }
    }
    if (dummyFlag != s_des)
    {
        return 0;
    }
}