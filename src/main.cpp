#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <optional>
#include <format>

#include "levenshtein_distance.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "symmetric_matrix.hpp"
#include "formulations.hpp"
#include "logger.hpp"
#include "ortools/base/init_google.h"
#include "ortools/init/init.h"

using namespace std;

void print_solution(optional<Solution> &solution, ofstream &output_file)
{
    if (solution.has_value())
    {
        string sanitized_string = solution.value().str.size() ? solution.value().str : "ε";

        output_file << "STATUS: " << solution.value().status << "\n"
                    << "SOLUTION: " << sanitized_string << "\n"
                    << "DISTANCE: " << solution.value().distance << "\n"
                    << "RUNTIME: " << solution.value().runtime << "s\n"
                    << "BEST LOWER BOUND: " << solution.value().lower_bound << "\n"
                    << "GAP: " << solution.value().gap << "\n";
    }
    else
    {
        output_file << "No solution was found.\n";
    }

    output_file << "\n";
    output_file.flush();
}

void print_if_inconsistent_solution(optional<Solution> &sol, vector<string> &input_strings, string_view label)
{
    if (!sol.has_value())
        return;

    double higher_distance = 0;

    for (string s : input_strings)
    {
        double distance = levenshtein_distance(sol.value().str, s);

        if (distance > higher_distance)
            higher_distance = distance;
    }

    if (higher_distance > sol.value().distance)
        Logger::error(format("\nINCONSISTENT SOLUTION - {}\n", label));
}

int main(int argc, char *argv[])
{
    string model = argv[3];

    ifstream input_file(argv[1]);
    ofstream output_file(argv[2], (model == "1" ? ios::out | ios::trunc : ios::out | ios::app));

    absl::SetStderrThreshold(absl::LogSeverityAtLeast::kError);
    InitGoogle(argv[0], &argc, &argv, true);

    string alphabet;
    input_file >> alphabet;
    StringCodec::set_alphabet(alphabet);

    int number_of_input_strings;
    input_file >> number_of_input_strings;

    vector<string> input_strings(number_of_input_strings);

    StringCodec codec;
    vector<vector<int>> encoded_input_strings(number_of_input_strings);

    int higher_str_size = numeric_limits<int>::min();
    int lower_str_size = numeric_limits<int>::max();

    for (int i = 0; i < number_of_input_strings; i++)
    {
        input_file >> input_strings[i];

        codec.encode(input_strings[i], encoded_input_strings[i]);

        int str_length = input_strings[i].size();

        if (str_length > higher_str_size)
            higher_str_size = str_length;

        if (str_length < lower_str_size)
            lower_str_size = str_length;
    }

    input_file.close();

    SymmetricMatrix<double> distance_matrix(input_strings.size());
    double diameter = numeric_limits<double>::lowest();

    for (int i = 0; i < input_strings.size(); i++)
    {
        for (int j = i; j < input_strings.size(); j++)
        {
            double distance = levenshtein_distance(
                encoded_input_strings[i],
                encoded_input_strings[j]);

            distance_matrix(i, j) = distance;

            if (distance > diameter)
                diameter = distance;
        }
    }

    double restricted_radius = numeric_limits<double>::max();
    for (int i = 0; i < input_strings.size(); i++)
    {
        double local_higher_distance = numeric_limits<double>::lowest();

        for (int j = 0; j < input_strings.size(); j++)
            if (distance_matrix(i, j) > local_higher_distance)
                local_higher_distance = distance_matrix(i, j);

        if (local_higher_distance < restricted_radius)
            restricted_radius = local_higher_distance;
    }

    const double SOL_UB = min(restricted_radius, CINS * higher_str_size);
    const double SOL_LB = diameter / 2;
    const int SOL_SIZE_UB = floor(lower_str_size + SOL_UB / CDEL);
    const int SOL_SIZE_LB = CSUB <= CINS ? lower_str_size : ceil(higher_str_size - SOL_UB / CINS);

    optional<Solution> solution;
    string model_name;

    if (model == "1")
    {
        model_name = "Hayashida and Koyano 2017 Gurobi";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::hayashida_and_koyano_2017(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, distance_matrix, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "2")
    {
        model_name = "MODIFIED Hayashida and Koyano 2017 Gurobi";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::improved_hayashida_and_koyano_2017(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, distance_matrix,
            higher_str_size, lower_str_size, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "3")
    {
        model_name = "Hayashida and Koyano 2017 CP-SAT";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::sat::hayashida_and_koyano_2017(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "4")
    {
        model_name = "MODIFIED Hayashida and Koyano 2017 CP-SAT";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::sat::improved_hayashida_and_koyano_2017(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, higher_str_size, lower_str_size, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }

    if (alphabet == "01" || alphabet == "10")
    {
        if (model == "5")
        {
            model_name = "Arbib et al. binary path Gurobi";
            Logger::info(format("Starting model {}...", model_name));
            output_file << format("({})", model_name) << endl;

            solution = operations_research::arbib_et_al_binary_path(
                SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
                SOL_UB, SOL_LB);

            print_solution(solution, output_file);
            Logger::success(format("Model {} execution completed.", model_name));
            // print_if_inconsistent_solution(solution, input_strings, model_name);
        }
        else if (model == "6")
        {
            model_name = "MODIFIED Arbib et al. binary path Gurobi";
            Logger::info(format("Starting model {}...", model_name));
            output_file << format("({})", model_name) << endl;

            solution = operations_research::improved_arbib_et_al_binary_path(
                SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
                higher_str_size, lower_str_size,
                SOL_UB, SOL_LB);

            print_solution(solution, output_file);
            Logger::success(format("Model {} execution completed.", model_name));
            // print_if_inconsistent_solution(solution, input_strings, model_name);
        }
        else if (model == "7")
        {
            model_name = "Arbib et al. binary path CP-SAT";
            Logger::info(format("Starting model {}...", model_name));
            output_file << format("({})", model_name) << endl;

            solution = operations_research::sat::arbib_et_al_binary_path(
                SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
                SOL_UB, SOL_LB);

            print_solution(solution, output_file);
            Logger::success(format("Model {} execution completed.", model_name));
            // print_if_inconsistent_solution(solution, input_strings, model_name);
        }
        else if (model == "8")
        {
            model_name = "MODIFIED Arbib et al. binary path CP-SAT";
            Logger::info(format("Starting model {}...", model_name));
            output_file << format("({})", model_name) << endl;

            solution = operations_research::sat::improved_arbib_et_al_binary_path(
                SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
                higher_str_size, lower_str_size,
                SOL_UB, SOL_LB);

            print_solution(solution, output_file);
            Logger::success(format("Model {} execution completed.", model_name));
            // print_if_inconsistent_solution(solution, input_strings, model_name);
        }
        else if (model == "9")
        {
            model_name = "Arbib et al. binary matching Gurobi";
            Logger::info(format("Starting model {}...", model_name));
            output_file << format("({})", model_name) << endl;

            solution = operations_research::arbib_et_al_binary_matching(
                SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
                SOL_UB, SOL_LB);

            print_solution(solution, output_file);
            Logger::success(format("Model {} execution completed.", model_name));
            // print_if_inconsistent_solution(solution, input_strings, model_name);
        }
        else if (model == "10")
        {
            model_name = "MODIFIED Arbib et al. binary matching Gurobi";
            Logger::info(format("Starting model {}...", model_name));
            output_file << format("({})", model_name) << endl;

            solution = operations_research::improved_arbib_et_al_binary_matching(
                SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
                higher_str_size, lower_str_size,
                SOL_UB, SOL_LB);

            print_solution(solution, output_file);
            Logger::success(format("Model {} execution completed.", model_name));
            // print_if_inconsistent_solution(solution, input_strings, model_name);
        }
        else if (model == "11")
        {
            model_name = "Arbib et al. binary matching CP-SAT";
            Logger::info(format("Starting model {}...", model_name));
            output_file << format("({})", model_name) << endl;

            solution = operations_research::sat::arbib_et_al_binary_matching(
                SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
                SOL_UB, SOL_LB);

            print_solution(solution, output_file);
            Logger::success(format("Model {} execution completed.", model_name));
            // print_if_inconsistent_solution(solution, input_strings, model_name);
        }
        else if (model == "12")
        {
            model_name = "MODIFIED Arbib et al. binary matching CP-SAT";
            Logger::info(format("Starting model {}...", model_name));
            output_file << format("({})", model_name) << endl;

            solution = operations_research::sat::improved_arbib_et_al_binary_matching(
                SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
                higher_str_size, lower_str_size,
                SOL_UB, SOL_LB);

            print_solution(solution, output_file);
            Logger::success(format("Model {} execution completed.", model_name));
            // print_if_inconsistent_solution(solution, input_strings, model_name);
        }
    }

    if (model == "13")
    {
        model_name = "Arbib et al. path Gurobi";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::arbib_et_al_path(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "14")
    {
        model_name = "MODIFIED Arbib et al. path Gurobi";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::improved_arbib_et_al_path(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, higher_str_size, lower_str_size, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "15")
    {
        model_name = "Arbib et al. matching Gurobi";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::arbib_et_al_matching(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "16")
    {
        model_name = "MODIFIED Arbib et al. matching Gurobi";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::improved_arbib_et_al_matching(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, higher_str_size, lower_str_size, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "17")
    {
        model_name = "Arbib et al. path CP-SAT";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::sat::arbib_et_al_path(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "18")
    {
        model_name = "MODIFIED Arbib et al. path CP-SAT";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::sat::improved_arbib_et_al_path(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, higher_str_size, lower_str_size, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "19")
    {
        model_name = "Arbib et al. matching CP-SAT";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::sat::arbib_et_al_matching(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }
    else if (model == "20")
    {
        model_name = "MODIFIED Arbib et al. matching CP-SAT";
        Logger::info(format("Starting model {}...", model_name));
        output_file << format("({})", model_name) << endl;

        solution = operations_research::sat::improved_arbib_et_al_matching(
            SOL_SIZE_UB, SOL_SIZE_LB, encoded_input_strings,
            alphabet, higher_str_size, lower_str_size, SOL_UB, SOL_LB);

        print_solution(solution, output_file);
        Logger::success(format("Model {} execution completed.", model_name));
        // print_if_inconsistent_solution(solution, input_strings, model_name);
    }

    output_file.close();

    return 0;
}
