#include <memory>

#include "formulations.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "symmetric_matrix.hpp"
#include "ortools/linear_solver/linear_solver.h"

using namespace std;

namespace operations_research
{
    optional<Solution> improved_hayashida_and_koyano_2017(
        int m,
        const vector<vector<int>> &s,
        const string_view &a,
        const SymmetricMatrix<double> &dm,
        int lower_str_size,
        int higher_str_size,
        int solution_lower_bound)
    {
        unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GUROBI"));
        if (!solver)
        {
            LOG(WARNING) << "Could not create solver GUROBI";
            return nullopt;
        }

        solver->set_time_limit(TIME_LIMIT_IN_MS);

        MPVariable *d = solver->MakeNumVar(solution_lower_bound, higher_str_size, "d");
        MPVariable *l = solver->MakeIntVar(lower_str_size, m, "l");

        vector<MPVariable *> t;
        solver->MakeIntVarArray(m, 0.0, a.size() - 1, "t", &t);

        vector<vector<vector<MPVariable *>>> x(s.size() + 1);
        vector<vector<vector<MPVariable *>>> y(s.size() + 1);
        vector<vector<vector<MPVariable *>>> z(s.size() + 1);
        vector<vector<vector<MPVariable *>>> g(s.size() + 1);
        vector<vector<vector<MPVariable *>>> h(s.size() + 1);

        for (int k = 0; k < s.size(); k++)
        {
            x[k].resize(s[k].size() + 1);
            y[k].resize(s[k].size() + 1);
            z[k].resize(s[k].size() + 1);
            g[k].resize(s[k].size() + 1);
            h[k].resize(s[k].size() + 1);

            for (int i = 0; i <= s[k].size(); i++)
            {
                x[k][i].resize(m + 1);
                y[k][i].resize(m + 1);
                z[k][i].resize(m + 1);
                g[k][i].resize(m + 1);
                h[k][i].resize(m + 1);

                for (int j = 0; j <= m; j++)
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

        LinearExpr diff = m - LinearExpr(l);

        for (int k = 0; k < s.size(); k++)
        {
            LinearExpr distance;

            size_t nk = s[k].size();

            for (int i = 1; i <= nk; i++)
            {
                distance += CDEL * LinearExpr(x[k][i][0]);

                for (int j = 1; j <= m; j++)
                {
                    distance += (CDEL * LinearExpr(x[k][i][j]) +
                                 CINS * LinearExpr(y[k][i][j]) +
                                 CSUB * LinearExpr(h[k][i][j]));
                }
            }

            for (int j = 1; j <= m; j++)
                distance += CINS * LinearExpr(y[k][0][j]);

            distance -= CINS * diff;

            solver->MakeRowConstraint(distance <= d);

            solver->MakeRowConstraint(
                LinearExpr(x[k][1][0]) +
                    LinearExpr(y[k][0][1]) +
                    LinearExpr(z[k][1][1]) ==
                1);

            solver->MakeRowConstraint(LinearExpr(x[k][nk][0]) == LinearExpr(y[k][nk][1]));

            solver->MakeRowConstraint(LinearExpr(x[k][1][m]) == LinearExpr(y[k][0][m]));

            solver->MakeRowConstraint(
                LinearExpr(x[k][nk][m]) +
                    LinearExpr(y[k][nk][m]) +
                    LinearExpr(z[k][nk][m]) ==
                1);

            for (int i = 1; i < nk; i++)
            {
                solver->MakeRowConstraint(
                    LinearExpr(x[k][i][0]) ==
                    LinearExpr(x[k][i + 1][0]) +
                        LinearExpr(y[k][i][1]) +
                        LinearExpr(z[k][i + 1][1]));

                solver->MakeRowConstraint(
                    LinearExpr(x[k][i][m]) +
                        LinearExpr(y[k][i][m]) +
                        LinearExpr(z[k][i][m]) ==
                    LinearExpr(x[k][i + 1][m]));

                for (int j = 1; j < m; j++)
                {
                    solver->MakeRowConstraint(
                        LinearExpr(x[k][i][j]) +
                            LinearExpr(y[k][i][j]) +
                            LinearExpr(z[k][i][j]) ==
                        LinearExpr(x[k][i + 1][j]) +
                            LinearExpr(y[k][i][j + 1]) +
                            LinearExpr(z[k][i + 1][j + 1]));
                }
            }

            for (int j = 1; j < m; j++)
            {
                solver->MakeRowConstraint(
                    LinearExpr(y[k][0][j]) ==
                    LinearExpr(x[k][1][j]) +
                        LinearExpr(y[k][0][j + 1]) +
                        LinearExpr(z[k][1][j + 1]));

                solver->MakeRowConstraint(
                    LinearExpr(x[k][nk][j]) +
                        LinearExpr(y[k][nk][j]) +
                        LinearExpr(z[k][nk][j]) ==
                    LinearExpr(y[k][nk][j + 1]));
            }

            for (int j = 1; j <= m; j++)
            {
                solver->MakeRowConstraint(
                    LinearExpr(y[k][nk][j]) >= 1.0 / m * (j - LinearExpr(l)));

                for (int i = 1; i <= nk; i++)
                {
                    solver->MakeRowConstraint(s[k][i - 1] - LinearExpr(t[j - 1]) <= a.size() * LinearExpr(g[k][i][j]));

                    solver->MakeRowConstraint(LinearExpr(t[j - 1]) - s[k][i - 1] <= a.size() * LinearExpr(g[k][i][j]));

                    solver->MakeRowConstraint(
                        LinearExpr(h[k][i][j]) >=
                        LinearExpr(z[k][i][j]) +
                            LinearExpr(g[k][i][j]) - 1);

                    solver->MakeRowConstraint(
                        LinearExpr(h[k][i][j]) <=
                        0.5 * (LinearExpr(z[k][i][j]) + LinearExpr(g[k][i][j])));
                }
            }
        }

        // triangle inequality
        for (int k1 = 0; k1 < s.size(); k1++)
        {
            LinearExpr distance_k1_t;
            size_t nk1 = s[k1].size();

            for (int i = 1; i <= nk1; i++)
            {
                distance_k1_t += CDEL * LinearExpr(x[k1][i][0]);

                for (int j = 1; j <= m; j++)
                {
                    distance_k1_t += CDEL * LinearExpr(x[k1][i][j]) +
                                     CINS * LinearExpr(y[k1][i][j]) +
                                     CSUB * LinearExpr(h[k1][i][j]);
                }
            }

            for (int j = 1; j <= m; j++)
                distance_k1_t += CINS * LinearExpr(y[k1][0][j]);

            distance_k1_t -= CINS * diff;

            for (int k2 = 0; k2 < s.size(); k2++)
            {
                if (k2 == k1)
                    continue;

                LinearExpr distance_k2_t;

                size_t nk2 = s[k2].size();

                for (int i = 1; i <= nk2; i++)
                {
                    distance_k2_t += CDEL * LinearExpr(x[k2][i][0]);

                    for (int j = 1; j <= m; j++)
                    {
                        distance_k2_t += CDEL * LinearExpr(x[k2][i][j]) +
                                         CINS * LinearExpr(y[k2][i][j]) +
                                         CSUB * LinearExpr(h[k2][i][j]);
                    }
                }

                for (int j = 1; j <= m; j++)
                    distance_k2_t += CINS * LinearExpr(y[k2][0][j]);

                distance_k2_t -= CINS * diff;

                solver->MakeRowConstraint(
                    dm(k1, k2) + distance_k1_t >= distance_k2_t);
            }

            for (int k2 = k1 + 1; k2 < s.size(); k2++)
            {
                LinearExpr distance_k2_t;

                size_t nk2 = s[k2].size();

                for (int i = 1; i <= nk2; i++)
                {
                    distance_k2_t += CDEL * LinearExpr(x[k2][i][0]);

                    for (int j = 1; j <= m; j++)
                    {
                        distance_k2_t += CDEL * LinearExpr(x[k2][i][j]) +
                                         CINS * LinearExpr(y[k2][i][j]) +
                                         CSUB * LinearExpr(h[k2][i][j]);
                    }
                }

                for (int j = 1; j <= m; j++)
                    distance_k2_t += CINS * LinearExpr(y[k2][0][j]);

                distance_k2_t -= CINS * diff;

                solver->MakeRowConstraint(
                    distance_k1_t + distance_k2_t >= dm(k1, k2));
            }
        }

        // additional constraints
        solver->MakeRowConstraint(higher_str_size - LinearExpr(d) <= LinearExpr(l));
        solver->MakeRowConstraint(LinearExpr(l) <= lower_str_size + LinearExpr(d));

        MPObjective *objective = solver->MutableObjective();
        objective->SetCoefficient(d, 1);
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
            d->solution_value(),
            status,
            absl::ToDoubleSeconds(solver->DurationSinceConstruction())};

        return solution;
    }
}