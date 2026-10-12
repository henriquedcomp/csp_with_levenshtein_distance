#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <optional>
#include <cstdlib>
#include <format>

#include "levenshtein_distance.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "symmetric_matrix.hpp"
#include "formulations.hpp"
#include "ortools/base/init_google.h"
#include "ortools/init/init.h"

using namespace std;

void print_solution(optional<Solution> &solution, ofstream &output_file, const string_view &title, bool last_line = false)
{
    output_file << format("({})", title) << "\n";

    if (solution.has_value())
    {
        string sanitized_string = solution.value().str.size() ? solution.value().str : "ε";

        output_file << "STATUS: " << solution.value().status << "\n"
                    << "SOLUTION: " << sanitized_string << "\n"
                    << "DISTANCE: " << solution.value().distance << "\n"
                    << "RUNTIME: " << solution.value().runtime << "s" << (!last_line ? "\n" : "");
    }
    else
    {
        output_file << "No solution was found." << (!last_line ? "\n" : "");
    }

    output_file << (!last_line ? "\n" : "");
}

int main(int argc, char *argv[])
{
    ifstream input_file(argv[1]);
    ofstream output_file(argv[2]);

    absl::SetStderrThreshold(absl::LogSeverityAtLeast::kInfo);
    InitGoogle(argv[0], &argc, &argv, true);

    string alphabet;
    int alphabet_size;

    input_file >> alphabet;
    alphabet_size = alphabet.size();

    StringCodec::set_alphabet(alphabet);

    int number_of_input_strings;

    input_file >> number_of_input_strings;

    vector<string> input_strings(number_of_input_strings);
    StringCodec codec;
    vector<vector<int>> encoded_input_strings(number_of_input_strings);

    int higher_str_size = numeric_limits<int>::min();
    int lower_str_size = numeric_limits<int>::max();
    int sum_str_size = 0;

    for (int i = 0; i < number_of_input_strings; i++)
    {
        input_file >> input_strings[i];

        codec.encode(input_strings[i], encoded_input_strings[i]);

        int str_length = input_strings[i].size();

        if (str_length > higher_str_size)
            higher_str_size = str_length;

        if (str_length < lower_str_size)
            lower_str_size = str_length;

        sum_str_size += str_length;
    }

    input_file.close();

    SymmetricMatrix<double> distance_matrix(input_strings.size());
    double diameter = numeric_limits<double>::min();
    int left_string_index, right_string_index;

    for (int i = 0; i < input_strings.size(); i++)
    {
        for (int j = i; j < input_strings.size(); j++)
        {
            double distance = levenshtein_distance(
                encoded_input_strings[i],
                encoded_input_strings[j]);

            distance_matrix(i, j) = distance;

            if (distance > diameter)
            {
                diameter = distance;
                left_string_index = i;
                right_string_index = j;
            }
        }
    }

    double lower_bound = diameter / 2;

    optional<Solution> solution = operations_research::hayashida_and_koyano_2017(
        floor(CSUB / CDEL * lower_str_size + higher_str_size),
        encoded_input_strings,
        alphabet,
        distance_matrix);

    print_solution(solution, output_file, "Hayashida and Koyano 2017 Gurobi");

    solution = operations_research::improved_hayashida_and_koyano_2017(
        floor(CSUB / CDEL * lower_str_size + higher_str_size),
        encoded_input_strings,
        alphabet,
        distance_matrix,
        lower_str_size,
        higher_str_size,
        lower_bound);

    print_solution(solution, output_file, "IMPROVED Hayashida and Koyano 2017 Gurobi");

    solution = operations_research::sat::hayashida_and_koyano_2017(
        sum_str_size,
        encoded_input_strings,
        alphabet,
        distance_matrix,
        higher_str_size,
        lower_str_size);

    print_solution(solution, output_file, "Hayashida and Koyano 2017 CP-SAT");

    solution = operations_research::sat::improved_hayashida_and_koyano_2017(
        floor(CSUB / CDEL * lower_str_size + higher_str_size),
        encoded_input_strings,
        alphabet,
        distance_matrix,
        lower_str_size,
        higher_str_size,
        lower_bound);

    print_solution(solution, output_file, "IMPROVED Hayashida and Koyano 2017 CP-SAT");

    if (alphabet == "01" || alphabet == "10")
    {
        solution = operations_research::arbib_et_al_binary_path(
            floor(CSUB / CDEL * lower_str_size + higher_str_size),
            encoded_input_strings,
            lower_str_size);

        print_solution(solution, output_file, "Arbib et al. binary path Gurobi");

        solution = operations_research::improved_arbib_et_al_binary_path(
            floor(CSUB / CDEL * lower_str_size + higher_str_size),
            encoded_input_strings,
            lower_str_size,
            distance_matrix,
            higher_str_size,
            lower_bound);

        print_solution(solution, output_file, "IMPROVED Arbib et al. binary path Gurobi");

        solution = operations_research::sat::arbib_et_al_binary_path(
            floor(CSUB / CDEL * lower_str_size + higher_str_size),
            encoded_input_strings,
            lower_str_size,
            higher_str_size);

        print_solution(solution, output_file, "Arbib et al. binary path CP-SAT");

        solution = operations_research::sat::improved_arbib_et_al_binary_path(
            floor(CSUB / CDEL * lower_str_size + higher_str_size),
            encoded_input_strings,
            lower_str_size,
            distance_matrix,
            higher_str_size,
            lower_bound);

        print_solution(solution, output_file, "IMPROVED Arbib et al. binary path CP-SAT");
    }

    solution = operations_research::arbib_et_al_path(
        floor(CSUB / CDEL * lower_str_size + higher_str_size),
        encoded_input_strings,
        alphabet,
        lower_str_size);

    print_solution(solution, output_file, "Arbib et al. path Gurobi");

    solution = operations_research::improved_arbib_et_al_path(
        floor(CSUB / CDEL * lower_str_size + higher_str_size),
        encoded_input_strings,
        alphabet,
        lower_str_size,
        distance_matrix,
        higher_str_size,
        lower_bound);

    print_solution(solution, output_file, "IMPROVED Arbib et al. path Gurobi");

    solution = operations_research::sat::arbib_et_al_path(
        floor(CSUB / CDEL * lower_str_size + higher_str_size),
        encoded_input_strings,
        alphabet,
        lower_str_size,
        higher_str_size);

    print_solution(solution, output_file, "Arbib et al. path CP-SAT");

    solution = operations_research::sat::improved_arbib_et_al_path(
        floor(CSUB / CDEL * lower_str_size + higher_str_size),
        encoded_input_strings,
        alphabet,
        lower_str_size,
        distance_matrix,
        higher_str_size,
        lower_bound);

    print_solution(solution, output_file, "IMPROVED Arbib et al. path CP-SAT", true);

    output_file.close();

    return 0;
}
