/*
    This script explores the tuples in C++, so they can be used poperly
    as dicts in Fitness Landscapes later.
*/

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

int is_in_dict(vector<char> desired, vector<tuple<vector<char>, int>> tup_vec)
{
    int s = tup_vec.size();
    int s_des = desired.size();

    for (int j = 0; j < s; j++)
    {
        auto seq = tup_vec[j];
        if (s_des != get<0>(seq).size())
        {
            continue; // Skip if sizes are different
        }

        bool all_equal = true; // Flag to check if all characters are equal

        for (int k = 0; k < s_des; k++)
        {
            if (desired[k] != get<0>(seq)[k])
            {
                all_equal = false;
                break; // Exit the loop if any character is not equal
            }
        }

        if (all_equal)
        {
            return j + 1;
        }
    }
    return 0;
}