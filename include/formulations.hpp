#pragma once

#include <string_view>
#include <optional>
#include <vector>

#include "symmetric_matrix.hpp"

using namespace std;

struct Solution
{
    string str;
    double distance;
    string status;
    double runtime;
};

namespace operations_research
{
    optional<Solution> hayashida_and_koyano_2017(
        int m,
        const vector<vector<int>> &s,
        const string_view &a,
        const SymmetricMatrix<double> &dm);

    namespace sat
    {
        optional<Solution> hayashida_and_koyano_2017(
            int m,
            const vector<vector<int>> &s,
            const string_view &a,
            const SymmetricMatrix<double> &dm);
    }
}