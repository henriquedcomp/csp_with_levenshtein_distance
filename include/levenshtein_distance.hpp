#pragma once

#include <vector>
#include <algorithm>

#include "config.hpp"

using namespace std;

double levenshtein_distance(const vector<int> &s, const vector<int> &t)
{
    vector<vector<double>> dp_matrix(2, vector<double>(t.size() + 1));

    for (int j = 0; j < t.size() + 1; j++)
    {
        dp_matrix[0][j] = j;
    }

    for (int i = 1; i < s.size() + 1; i++)
    {
        dp_matrix[1][0] = i;

        for (int j = 1; j < t.size() + 1; j++)
        {
            dp_matrix[1][j] = min({dp_matrix[0][j - 1] + (s[i - 1] == t[j - 1] ? 0 : CSUB),
                                   dp_matrix[0][j] + CINS,
                                   dp_matrix[1][j - 1] + CDEL});
        }

        swap(dp_matrix[0], dp_matrix[1]);
    }

    return dp_matrix[0][t.size()];
};

// I'll keep it here for possible future use (backtracking)
// double levenshtein_distance(const vector<int> &s, const vector<int> &t)
// {
//     vector<vector<double>> dp_matrix(s.size() + 1, vector<double>(t.size() + 1));

//     for (int i = 0; i < s.size() + 1; i++)
//     {
//         dp_matrix[i][0] = i;
//     }

//     for (int j = 1; j < t.size() + 1; j++)
//     {
//         dp_matrix[0][j] = j;
//     }

//     for (int i = 1; i < s.size() + 1; i++)
//     {
//         for (int j = 1; j < t.size() + 1; j++)
//         {
//             dp_matrix[i][j] = min({dp_matrix[i - 1][j - 1] + SUBS_MATRIX[s[i - 1]][t[j - 1]],
//                                    dp_matrix[i - 1][j] + INDEL_COST,
//                                    dp_matrix[i][j - 1] + INDEL_COST});
//         }
//     }

//     return dp_matrix[s.size()][t.size()];
// };
