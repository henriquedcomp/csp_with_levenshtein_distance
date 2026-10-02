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
        int m,                            // upper bound on solution size
        const vector<vector<int>> &s,     // encoded input strings
        const string_view &a,             // alphabet
        const SymmetricMatrix<double> &dm // distance matrix between input strings
    );

    optional<Solution> improved_hayashida_and_koyano_2017(
        int m,
        const vector<vector<int>> &s,
        const string_view &a,
        const SymmetricMatrix<double> &dm,
        int lower_str_size,
        int higher_str_size,
        int solution_lower_bound);

    optional<Solution> arbib_et_al_binary_path(
        int n,                                 // upper bound on solution size
        const vector<vector<int>> &stringpool, // input strings in 0-1 alphabet
        int p                                  // lower str size
    );

    optional<Solution> improved_arbib_et_al_binary_path(
        int n,
        const vector<vector<int>> &stringpool,
        int p,
        const SymmetricMatrix<double> &dm,
        int higher_str_size,
        int solution_lower_bound);

    optional<Solution> arbib_et_al_path(
        int n,
        const vector<vector<int>> &stringpool,
        const string_view &a,
        int p);

    optional<Solution> improved_arbib_et_al_path(
        int n,
        const vector<vector<int>> &stringpool,
        const string_view &a,
        int p,
        const SymmetricMatrix<double> &dm,
        int higher_str_size,
        int solution_lower_bound);

    namespace sat
    {
        optional<Solution> hayashida_and_koyano_2017(
            int m,
            const vector<vector<int>> &s,
            const string_view &a,
            const SymmetricMatrix<double> &dm);

        optional<Solution> improved_hayashida_and_koyano_2017(
            int m,
            const vector<vector<int>> &s,
            const string_view &a,
            const SymmetricMatrix<double> &dm,
            int lower_str_size,
            int higher_str_size,
            int solution_lower_bound);

        optional<Solution> arbib_et_al_binary_path(
            int n,
            const vector<vector<int>> &stringpool,
            int p);

        optional<Solution> improved_arbib_et_al_binary_path(
            int n,
            const vector<vector<int>> &stringpool,
            int p,
            const SymmetricMatrix<double> &dm,
            int higher_str_size,
            int solution_lower_bound);

        optional<Solution> arbib_et_al_path(
            int n,
            const vector<vector<int>> &stringpool,
            const string_view &a,
            int p);

        optional<Solution> improved_arbib_et_al_path(
            int n,
            const vector<vector<int>> &stringpool,
            const string_view &a,
            int p,
            const SymmetricMatrix<double> &dm,
            int higher_str_size,
            int solution_lower_bound);
    }
}