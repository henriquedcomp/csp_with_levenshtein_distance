#include <memory>

#include "formulations.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "symmetric_matrix.hpp"
#include "ortools/linear_solver/linear_solver.h"

using namespace std;

namespace operations_research
{
    optional<Solution> hayashida_and_koyano_2017(
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
        const SymmetricMatrix<double> &distance_matrix,
        double sol_ub,
        double sol_lb)
    {
        unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GUROBI"));
        if (!solver)
        {
            LOG(WARNING) << "Could not create solver GUROBI";
            return nullopt;
        }

        solver->set_time_limit(TIME_LIMIT_IN_MS);

        MPVariable *r = solver->MakeNumVar(sol_lb, sol_ub, "r");
        MPVariable *l = solver->MakeIntVar(sol_size_lb, sol_size_ub, "l");

        vector<MPVariable *> t;
        solver->MakeIntVarArray(sol_size_ub, 0.0, alphabet.size() - 1, "t", &t);

        vector<vector<vector<MPVariable *>>> x(input_strings.size() + 1);
        vector<vector<vector<MPVariable *>>> y(input_strings.size() + 1);
        vector<vector<vector<MPVariable *>>> z(input_strings.size() + 1);
        vector<vector<vector<MPVariable *>>> g(input_strings.size() + 1);
        vector<vector<vector<MPVariable *>>> h(input_strings.size() + 1);

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
                        string x_var_name = format("x_{}_{}", k, i);
                        x[k][i][j] = solver->MakeBoolVar(x_var_name);

                        if (j > 0)
                        {
                            string z_var_name = format("z_{}_{}", k, i);
                            string g_var_name = format("g_{}_{}", k, i);
                            string h_var_name = format("h_{}_{}", k, i);

                            z[k][i][j] = solver->MakeBoolVar(z_var_name);
                            g[k][i][j] = solver->MakeBoolVar(g_var_name);
                            h[k][i][j] = solver->MakeBoolVar(h_var_name);
                        }
                    }

                    if (j > 0)
                    {
                        string y_var_name = format("y_{}_{}", k, i);
                        y[k][i][j] = solver->MakeBoolVar(y_var_name);
                    }
                }
            }
        }

        LinearExpr diff = sol_size_ub - LinearExpr(l);

        for (int k = 0; k < input_strings.size(); k++)
        {
            LinearExpr distance;

            size_t nk = input_strings[k].size();

            for (int i = 1; i <= nk; i++)
            {
                distance += CDEL * LinearExpr(x[k][i][0]);

                for (int j = 1; j <= sol_size_ub; j++)
                {
                    distance += (CDEL * LinearExpr(x[k][i][j]) +
                                 CINS * LinearExpr(y[k][i][j]) +
                                 CSUB * LinearExpr(h[k][i][j]));
                }
            }

            for (int j = 1; j <= sol_size_ub; j++)
                distance += CINS * LinearExpr(y[k][0][j]);

            distance -= CINS * diff;

            solver->MakeRowConstraint(distance <= r); // (4.2)

            solver->MakeRowConstraint(
                LinearExpr(x[k][1][0]) +
                    LinearExpr(y[k][0][1]) +
                    LinearExpr(z[k][1][1]) ==
                1); // (4.3)

            solver->MakeRowConstraint(LinearExpr(x[k][nk][0]) == LinearExpr(y[k][nk][1])); // (4.5)

            solver->MakeRowConstraint(LinearExpr(x[k][1][sol_size_ub]) == LinearExpr(y[k][0][sol_size_ub])); // (4.7)

            solver->MakeRowConstraint(
                LinearExpr(x[k][nk][sol_size_ub]) +
                    LinearExpr(y[k][nk][sol_size_ub]) +
                    LinearExpr(z[k][nk][sol_size_ub]) ==
                1); // (4.11)

            for (int i = 1; i < nk; i++)
            {
                solver->MakeRowConstraint(
                    LinearExpr(x[k][i][0]) ==
                    LinearExpr(x[k][i + 1][0]) +
                        LinearExpr(y[k][i][1]) +
                        LinearExpr(z[k][i + 1][1])); // (4.4)

                solver->MakeRowConstraint(
                    LinearExpr(x[k][i][sol_size_ub]) +
                        LinearExpr(y[k][i][sol_size_ub]) +
                        LinearExpr(z[k][i][sol_size_ub]) ==
                    LinearExpr(x[k][i + 1][sol_size_ub])); // (4.10)

                for (int j = 1; j < sol_size_ub; j++)
                {
                    solver->MakeRowConstraint(
                        LinearExpr(x[k][i][j]) +
                            LinearExpr(y[k][i][j]) +
                            LinearExpr(z[k][i][j]) ==
                        LinearExpr(x[k][i + 1][j]) +
                            LinearExpr(y[k][i][j + 1]) +
                            LinearExpr(z[k][i + 1][j + 1])); // (4.8)
                }
            }

            for (int j = 1; j < sol_size_ub; j++)
            {
                solver->MakeRowConstraint(
                    LinearExpr(y[k][0][j]) ==
                    LinearExpr(x[k][1][j]) +
                        LinearExpr(y[k][0][j + 1]) +
                        LinearExpr(z[k][1][j + 1])); // (4.6)

                solver->MakeRowConstraint(
                    LinearExpr(x[k][nk][j]) +
                        LinearExpr(y[k][nk][j]) +
                        LinearExpr(z[k][nk][j]) ==
                    LinearExpr(y[k][nk][j + 1])); // (4.9)
            }

            for (int j = 1; j <= sol_size_ub; j++)
            {
                solver->MakeRowConstraint(
                    LinearExpr(y[k][nk][j]) >= 1.0 / sol_size_ub * (j - LinearExpr(l))); // (4.12)

                for (int i = 1; i <= nk; i++)
                {
                    solver->MakeRowConstraint(input_strings[k][i - 1] - LinearExpr(t[j - 1]) <=
                                              alphabet.size() * LinearExpr(g[k][i][j])); // (4.13)

                    solver->MakeRowConstraint(LinearExpr(t[j - 1]) - input_strings[k][i - 1] <=
                                              alphabet.size() * LinearExpr(g[k][i][j])); // (4.14)

                    solver->MakeRowConstraint(
                        LinearExpr(h[k][i][j]) >=
                        LinearExpr(z[k][i][j]) +
                            LinearExpr(g[k][i][j]) - 1); // (4.15)

                    solver->MakeRowConstraint(
                        LinearExpr(h[k][i][j]) <=
                        0.5 * (LinearExpr(z[k][i][j]) + LinearExpr(g[k][i][j]))); // (4.16)
                }
            }
        }

        // triangle inequality
        for (int k1 = 0; k1 < input_strings.size(); k1++)
        {
            LinearExpr distance_k1_t;
            size_t nk1 = input_strings[k1].size();

            for (int i = 1; i <= nk1; i++)
            {
                distance_k1_t += CDEL * LinearExpr(x[k1][i][0]);

                for (int j = 1; j <= sol_size_ub; j++)
                {
                    distance_k1_t += CDEL * LinearExpr(x[k1][i][j]) +
                                     CINS * LinearExpr(y[k1][i][j]) +
                                     CSUB * LinearExpr(h[k1][i][j]);
                }
            }

            for (int j = 1; j <= sol_size_ub; j++)
                distance_k1_t += CINS * LinearExpr(y[k1][0][j]);

            distance_k1_t -= CINS * diff;

            for (int k2 = 0; k2 < input_strings.size(); k2++)
            {
                if (k2 == k1)
                    continue;

                LinearExpr distance_k2_t;

                size_t nk2 = input_strings[k2].size();

                for (int i = 1; i <= nk2; i++)
                {
                    distance_k2_t += CDEL * LinearExpr(x[k2][i][0]);

                    for (int j = 1; j <= sol_size_ub; j++)
                    {
                        distance_k2_t += CDEL * LinearExpr(x[k2][i][j]) +
                                         CINS * LinearExpr(y[k2][i][j]) +
                                         CSUB * LinearExpr(h[k2][i][j]);
                    }
                }

                for (int j = 1; j <= sol_size_ub; j++)
                    distance_k2_t += CINS * LinearExpr(y[k2][0][j]);

                distance_k2_t -= CINS * diff;

                solver->MakeRowConstraint(
                    distance_matrix(k1, k2) + distance_k1_t >= distance_k2_t); // (4.17)
            }

            for (int k2 = k1 + 1; k2 < input_strings.size(); k2++)
            {
                LinearExpr distance_k2_t;

                size_t nk2 = input_strings[k2].size();

                for (int i = 1; i <= nk2; i++)
                {
                    distance_k2_t += CDEL * LinearExpr(x[k2][i][0]);

                    for (int j = 1; j <= sol_size_ub; j++)
                    {
                        distance_k2_t += CDEL * LinearExpr(x[k2][i][j]) +
                                         CINS * LinearExpr(y[k2][i][j]) +
                                         CSUB * LinearExpr(h[k2][i][j]);
                    }
                }

                for (int j = 1; j <= sol_size_ub; j++)
                    distance_k2_t += CINS * LinearExpr(y[k2][0][j]);

                distance_k2_t -= CINS * diff;

                solver->MakeRowConstraint(
                    distance_k1_t + distance_k2_t >= distance_matrix(k1, k2)); // (4.18)
            }
        }

        MPObjective *objective = solver->MutableObjective();
        objective->SetCoefficient(r, 1);
        objective->SetMinimization();

        const MPSolver::ResultStatus result_status = solver->Solve();
        string status = "OPTIMAL";

        if (result_status != MPSolver::OPTIMAL)
        {
            LOG(INFO) << "The solver could not find an optimal solution in specified duration!";
            if (result_status == MPSolver::FEASIBLE)
            {
                status = "FEASIBLE";
                LOG(INFO) << "A potentially suboptimal solution was found";
            }
            else
            {
                status = "";
                LOG(INFO) << "The solver could not solve the problem.";
                LOG(WARNING) << "Time elapsed: " << absl::ToDoubleSeconds(solver->DurationSinceConstruction()) << "s";
                return nullopt;
            }
        }

        vector<int> encoded_solution_str;
        string decoded_solution_str;

        encoded_solution_str.reserve(t.size());

        int i = 0;
        for (MPVariable *var : t)
        {
            if (i++ >= l->solution_value())
                break;

            int value = static_cast<int>(var->solution_value());

            encoded_solution_str.push_back(value);
        }

        StringCodec codec;
        codec.decode(encoded_solution_str, decoded_solution_str);

        Solution solution = {
            decoded_solution_str,
            r->solution_value(),
            status,
            absl::ToDoubleSeconds(solver->DurationSinceConstruction()),
            objective->BestBound(),
            (objective->Value() - objective->BestBound()) / objective->Value() * 100.0};

        return solution;
    }
}