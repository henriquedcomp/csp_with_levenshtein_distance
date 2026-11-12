#include "formulations.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "symmetric_matrix.hpp"
#include "ortools/sat/cp_model.h"

using namespace std;

namespace operations_research
{
    namespace sat
    {
        optional<Solution> improved_hayashida_and_koyano_2017(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            const string_view &alphabet,
            int higher_str_size,
            int lower_str_size,
            double sol_ub,
            double sol_lb)
        {
            CpModelBuilder cp_model;

            const Domain solution_domain(sol_lb, sol_ub);
            const IntVar r = cp_model.NewIntVar(solution_domain);

            const Domain solution_length_domain(sol_size_lb, sol_size_ub);
            const IntVar l = cp_model.NewIntVar(solution_length_domain);

            const Domain alphabet_domain(0, alphabet.size() - 1);
            vector<IntVar> t(sol_size_ub);

            for (int i = 0; i < sol_size_ub; i++)
                t[i] = cp_model.NewIntVar(alphabet_domain);

            vector<vector<vector<BoolVar>>> x(input_strings.size() + 1);
            vector<vector<vector<BoolVar>>> y(input_strings.size() + 1);
            vector<vector<vector<BoolVar>>> z(input_strings.size() + 1);
            vector<vector<vector<BoolVar>>> g(input_strings.size() + 1);
            vector<vector<vector<BoolVar>>> h(input_strings.size() + 1);

            for (int k = 0; k < input_strings.size(); k++)
            {
                x[k].resize(input_strings[k].size() + 1);
                y[k].resize(input_strings[k].size() + 1);
                z[k].resize(input_strings[k].size() + 1);
                g[k].resize(input_strings[k].size() + 1);
                h[k].resize(input_strings[k].size() + 1);

                for (int i = 0; i <= input_strings[k].size(); i++)
                {
                    x[k][i].resize(sol_size_ub + 1);
                    y[k][i].resize(sol_size_ub + 1);
                    z[k][i].resize(sol_size_ub + 1);
                    g[k][i].resize(sol_size_ub + 1);
                    h[k][i].resize(sol_size_ub + 1);

                    for (int j = 0; j <= sol_size_ub; j++)
                    {
                        if (i > 0)
                        {
                            x[k][i][j] = cp_model.NewBoolVar();

                            if (j > 0)
                            {
                                z[k][i][j] = cp_model.NewBoolVar();
                                g[k][i][j] = cp_model.NewBoolVar();
                                h[k][i][j] = cp_model.NewBoolVar();
                            }
                        }

                        if (j > 0)
                            y[k][i][j] = cp_model.NewBoolVar();
                    }
                }
            }

            LinearExpr diff = sol_size_ub - l;

            for (int k = 0; k < input_strings.size(); k++)
            {
                LinearExpr distance;

                size_t nk = input_strings[k].size();

                for (int i = 1; i <= nk; i++)
                {
                    distance += CDEL * x[k][i][0];

                    for (int j = 1; j <= sol_size_ub; j++)
                    {
                        distance += (CDEL * x[k][i][j] +
                                     CINS * y[k][i][j] +
                                     CSUB * h[k][i][j]);
                    }
                }

                for (int j = 1; j <= sol_size_ub; j++)
                    distance += CINS * y[k][0][j];

                distance -= CINS * diff;

                cp_model.AddLessOrEqual(distance, r); // (4.2)

                cp_model.AddExactlyOne({x[k][1][0], y[k][0][1], z[k][1][1]});                                  // (4.3)
                cp_model.AddEquality(x[k][nk][0], y[k][nk][1]);                                                // (4.5)
                cp_model.AddEquality(x[k][1][sol_size_ub], y[k][0][sol_size_ub]);                              // (4.7)
                cp_model.AddExactlyOne({x[k][nk][sol_size_ub], y[k][nk][sol_size_ub], z[k][nk][sol_size_ub]}); // (4.11)

                for (int i = 1; i < nk; i++)
                {
                    cp_model.AddEquality(x[k][i][0], x[k][i + 1][0] + y[k][i][1] + z[k][i + 1][1]); // (4.4)

                    cp_model.AddEquality(x[k][i][sol_size_ub] + y[k][i][sol_size_ub] + z[k][i][sol_size_ub], x[k][i + 1][sol_size_ub]); // (4.10)

                    for (int j = 1; j < sol_size_ub; j++)
                    {
                        cp_model.AddEquality(
                            x[k][i][j] + y[k][i][j] + z[k][i][j],
                            x[k][i + 1][j] + y[k][i][j + 1] + z[k][i + 1][j + 1]); // (4.8)
                    }
                }

                for (int j = 1; j < sol_size_ub; j++)
                {
                    cp_model.AddEquality(y[k][0][j], x[k][1][j] + y[k][0][j + 1] + z[k][1][j + 1]); // (4.6)
                    cp_model.AddEquality(x[k][nk][j] + y[k][nk][j] + z[k][nk][j], y[k][nk][j + 1]); // (4.9)
                }

                for (int j = 1; j <= sol_size_ub; j++)
                {
                    cp_model.AddLessOrEqual(j, l).OnlyEnforceIf(Not(y[k][nk][j])); // (4.26)

                    for (int i = 1; i <= nk; i++)
                    {
                        cp_model.AddEquality(input_strings[k][i - 1], t[j - 1]).OnlyEnforceIf(Not(g[k][i][j])); // (4.27)
                        cp_model.AddEquality(h[k][i][j], 1).OnlyEnforceIf({z[k][i][j], g[k][i][j]});            // (4.28)
                    }
                }
            }

            // additional constraints
            cp_model.AddLessOrEqual(CINS * higher_str_size - r, l * CINS); // (4.24)
            cp_model.AddLessOrEqual(l * CDEL, lower_str_size * CDEL + r);  // (4.25)

            cp_model.Minimize(r);

            Model model;
            SatParameters parameters;

            parameters.set_max_time_in_seconds(TIME_LIMIT_IN_S);
            model.Add(NewSatParameters(parameters));

            const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
            LOG(INFO) << CpSolverResponseStats(response);

            string status = "OPTIMAL";

            if (response.status() != CpSolverStatus::OPTIMAL)
            {
                LOG(INFO) << "The solver could not find an optimal solution in specified duration!";

                if (response.status() == CpSolverStatus::FEASIBLE)
                {
                    status = "FEASIBLE";
                    LOG(INFO) << "A potentially suboptimal solution was found";
                }
                else
                {
                    status = "";
                    LOG(INFO) << "The solver could not solve the problem.";
                    return nullopt;
                }
            }

            vector<int> encoded_solution_str;
            string decoded_solution_str;

            encoded_solution_str.reserve(t.size());

            int i = 0;
            for (IntVar var : t)
            {
                if (i++ >= SolutionIntegerValue(response, l))
                    break;

                int value = SolutionIntegerValue(response, var);

                encoded_solution_str.push_back(value);
            }

            StringCodec codec;
            codec.decode(encoded_solution_str, decoded_solution_str);

            Solution solution = {
                decoded_solution_str,
                static_cast<double>(SolutionIntegerValue(response, r)),
                status,
                response.wall_time(),
                response.best_objective_bound(),
                (response.objective_value() - response.best_objective_bound()) / response.objective_value() * 100.0};

            return solution;
        }
    }
}