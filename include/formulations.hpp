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
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
        const SymmetricMatrix<double> &distance_matrix,
        double sol_ub,
        double sol_lb);

    optional<Solution> improved_hayashida_and_koyano_2017(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
        const SymmetricMatrix<double> &distance_matrix,
        int higher_str_size,
        int lower_str_size,
        double sol_ub,
        double sol_lb);

    optional<Solution> arbib_et_al_binary_path(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        double sol_ub,
        double sol_lb);

    optional<Solution> arbib_et_al_binary_matching(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        double sol_ub,
        double sol_lb);

    optional<Solution> improved_arbib_et_al_binary_path(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        int higher_str_size,
        int lower_str_size,
        double sol_ub,
        double sol_lb);

    optional<Solution> improved_arbib_et_al_binary_matching(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        int higher_str_size,
        int lower_str_size,
        double sol_ub,
        double sol_lb);

    optional<Solution> arbib_et_al_path(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
        double sol_ub,
        double sol_lb);

    optional<Solution> arbib_et_al_matching(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
        double sol_ub,
        double sol_lb);

    optional<Solution> improved_arbib_et_al_path(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
        int higher_str_size,
        int lower_str_size,
        double sol_ub,
        double sol_lb);

    optional<Solution> improved_arbib_et_al_matching(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
        int higher_str_size,
        int lower_str_size,
        double sol_ub,
        double sol_lb);

    namespace sat
    {
        optional<Solution> hayashida_and_koyano_2017(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            const string_view &alphabet,
            double sol_ub,
            double sol_lb);

        optional<Solution> improved_hayashida_and_koyano_2017(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            const string_view &alphabet,
            int higher_str_size,
            int lower_str_size,
            double sol_ub,
            double sol_lb);

        optional<Solution> arbib_et_al_binary_path(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            double sol_ub,
            double sol_lb);

        optional<Solution> arbib_et_al_binary_matching(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            double sol_ub,
            double sol_lb);

        optional<Solution> improved_arbib_et_al_binary_path(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            int higher_str_size,
            int lower_str_size,
            double sol_ub,
            double sol_lb);

        optional<Solution> improved_arbib_et_al_binary_matching(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            int higher_str_size,
            int lower_str_size,
            double sol_ub,
            double sol_lb);

        optional<Solution> arbib_et_al_path(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            const string_view &alphabet,
            double sol_ub,
            double sol_lb);

        optional<Solution> arbib_et_al_matching(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            const string_view &alphabet,
            double sol_ub,
            double sol_lb);

        optional<Solution> improved_arbib_et_al_path(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            const string_view &alphabet,
            int higher_str_size,
            int lower_str_size,
            double sol_ub,
            double sol_lb);

        optional<Solution> improved_arbib_et_al_matching(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            const string_view &alphabet,
            int higher_str_size,
            int lower_str_size,
            double sol_ub,
            double sol_lb);
    }
}