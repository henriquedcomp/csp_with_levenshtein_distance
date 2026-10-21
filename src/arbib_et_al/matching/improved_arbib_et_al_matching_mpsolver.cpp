#include <memory>
#include <map>

#include "formulations.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "offset_vector.hpp"
#include "graph_construction.hpp"

#include "ortools/linear_solver/linear_solver.h"

using namespace std;

using StringArc = tuple<int, int, int>;
using StringIndelArc = tuple<int, int>;

namespace operations_research
{
    optional<Solution> improved_arbib_et_al_matching(
        int n,
        const vector<vector<int>> &stringpool,
        const string_view &a,
        int p,
        const SymmetricMatrix<double> &dm,
        int higher_str_size,
        double solution_lower_bound)
    {
        unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GUROBI"));
        if (!solver)
        {
            LOG(WARNING) << "Could not create solver GUROBI";
            return nullopt;
        }

        solver->set_time_limit(TIME_LIMIT_IN_MS);

        const double solution_ub = CSUB * p + CDEL * (higher_str_size - p);
        MPVariable *d = solver->MakeNumVar(solution_lower_bound, solution_ub, "d");

        vector<MPVariable *> x;
        solver->MakeNumVarArray(n, 0, a.size() - 1, "x", &x);

        OffsetVector<MPVariable *> z(1, n);
        for (int j = 1; j <= n; j++)
        {
            string variable_name = format("z_{}", j);
            z[j] = solver->MakeBoolVar(variable_name);
        }

        for (int j = 1; j < n; j++)
            solver->MakeRowConstraint(LinearExpr(z[j + 1]) <= LinearExpr(z[j])); // (18)

        solver->MakeRowConstraint(LinearExpr(z[1]) == 1.0);

        map<StringArc, MPVariable *> yF; // free substitutions
        map<StringArc, MPVariable *> yC; // costly substitutions

        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int h = 1; h <= stringpool[k].size(); h++)
            {
                for (int j = 1; j <= n; j++)
                {
                    string free_variable_name = format("yF_{}_{}_{}", k, h, j);
                    string costly_variable_name = format("yC_{}_{}_{}", k, h, j);

                    yF[StringArc(k, h, j)] = solver->MakeBoolVar(free_variable_name);
                    yC[StringArc(k, h, j)] = solver->MakeBoolVar(costly_variable_name);
                }
            }
        }

        map<StringIndelArc, MPVariable *> yDel;
        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int h = 1; h <= stringpool[k].size(); h++)
            {
                string variable_name = format("yDel_{}_{}", k, h);

                yDel[StringIndelArc(k, h)] = solver->MakeBoolVar(variable_name);
            }
        }

        map<StringIndelArc, MPVariable *> yIns;
        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int j = 1; j <= n; j++)
            {
                string variable_name = format("yIns_{}_{}", k, j);

                yIns[StringIndelArc(k, j)] = solver->MakeBoolVar(variable_name);
            }
        }

        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int h = 1; h <= stringpool[k].size(); h++)
            {
                LinearExpr sum = LinearExpr(yDel[StringIndelArc(k, h)]);

                for (int j = 1; j <= n; j++)
                    sum += LinearExpr(yF[StringArc(k, h, j)]) +
                           LinearExpr(yC[StringArc(k, h, j)]);

                solver->MakeRowConstraint(sum == 1.0);
            }
        }

        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int h = 1; h <= stringpool[k].size(); h++)
            {
                LinearExpr sumJ, sumH, sumL;

                for (int j = 1; j <= n; j++)
                    sumJ += LinearExpr(yF[StringArc(k, h, j)]) +
                            LinearExpr(yC[StringArc(k, h, j)]);

                for (int u = h + 1; u <= stringpool[k].size(); u++)
                    sumH += LinearExpr(yF[StringArc(k, u, 1)]) +
                            LinearExpr(yC[StringArc(k, u, 1)]);

                for (int u = 1; u < h; u++)
                    sumL += LinearExpr(yF[StringArc(k, u, n)]) +
                            LinearExpr(yC[StringArc(k, u, n)]);

                solver->MakeRowConstraint(sumJ + sumH + sumL <= 1.0);
            }
        }

        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int j = 1; j <= n; j++)
            {
                LinearExpr sum = LinearExpr(yIns[StringIndelArc(k, j)]);

                for (int h = 1; h <= stringpool[k].size(); h++)
                    sum += LinearExpr(yF[StringArc(k, h, j)]) +
                           LinearExpr(yC[StringArc(k, h, j)]);

                solver->MakeRowConstraint(sum == z[j]);
            }
        }

        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int j = 1; j <= n; j++)
            {
                LinearExpr sum_h, sumH, sumL;

                for (int h = 1; h <= stringpool[k].size(); h++)
                    sum_h += LinearExpr(yF[StringArc(k, h, j)]) +
                             LinearExpr(yC[StringArc(k, h, j)]);

                for (int v = j + 1; v <= n; v++)
                    sumH += LinearExpr(yF[StringArc(k, 1, v)]) +
                            LinearExpr(yC[StringArc(k, 1, v)]);

                for (int v = 1; v < j; v++)
                    sumL += LinearExpr(yF[StringArc(k, stringpool[k].size(), v)]) +
                            LinearExpr(yC[StringArc(k, stringpool[k].size(), v)]);

                solver->MakeRowConstraint(sum_h + sumH + sumL <= 1.0);
            }
        }

        // nocross constraint reimagined
        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int u = 1; u <= stringpool[k].size(); u++)
            {
                for (int v = 1; v <= n; v++)
                {
                    for (int j = v; j <= n; j++)
                    {
                        LinearExpr sum = LinearExpr(yF[StringArc(k, u, v)]) +
                                         LinearExpr(yC[StringArc(k, u, v)]);

                        for (int h = 1; h < u; h++)
                        {
                            sum += LinearExpr(yF[StringArc(k, h, j)]) +
                                   LinearExpr(yC[StringArc(k, h, j)]);
                        }

                        solver->MakeRowConstraint(sum <= 1.0);
                    }

                    for (int h = u; h <= stringpool[k].size(); h++)
                    {
                        LinearExpr sum = LinearExpr(yF[StringArc(k, u, v)]) +
                                         LinearExpr(yC[StringArc(k, u, v)]);

                        for (int j = 1; j < v; j++)
                        {
                            sum += LinearExpr(yF[StringArc(k, h, j)]) +
                                   LinearExpr(yC[StringArc(k, h, j)]);
                        }

                        solver->MakeRowConstraint(sum <= 1.0);
                    }
                }
            }
        }

        for (int k = 0; k < stringpool.size(); k++)
        {
            for (int h = 1; h <= stringpool[k].size(); h++)
            {
                int lam = stringpool[k][h - 1];

                for (int j = 1; j <= n; j++)
                {
                    string name = format("e_{}_{}_{}", k, h, j);
                    MPVariable *e = solver->MakeBoolVar(name);

                    solver->MakeRowConstraint(lam - LinearExpr(x[j - 1]) <= a.size() * LinearExpr(e));

                    solver->MakeRowConstraint(LinearExpr(x[j - 1]) - lam <= a.size() * LinearExpr(e));

                    solver->MakeRowConstraint(LinearExpr(yF[StringArc(k, h, j)]) + LinearExpr(e) <= 1.0);
                }
            }
        }

        for (int k = 0; k < stringpool.size(); k++)
        {
            LinearExpr distance;

            for (int h = 1; h <= stringpool[k].size(); h++)
            {
                distance += CDEL * LinearExpr(yDel[StringIndelArc(k, h)]);

                for (int j = 1; j <= n; j++)
                    distance += CSUB * LinearExpr(yC[StringArc(k, h, j)]);
            }

            for (int j = 1; j <= n; j++)
                distance += CINS * LinearExpr(yIns[StringIndelArc(k, j)]);

            solver->MakeRowConstraint(distance <= d);
        }

        // triangle inequality - it makes results even worse
        // for (int i = 0; i < stringpool.size(); i++)
        // {
        //     LinearExpr distance_i_x;

        //     for (int h = 1; h <= stringpool[i].size(); h++)
        //     {
        //         distance_i_x += CDEL * LinearExpr(yDel[StringIndelArc(i, h)]);

        //         for (int j = 1; j <= n; j++)
        //             distance_i_x += CSUB * LinearExpr(yC[StringArc(i, h, j)]);
        //     }

        //     for (int j = 1; j <= n; j++)
        //         distance_i_x += CINS * LinearExpr(yIns[StringIndelArc(i, j)]);

        //     for (int j = 0; j < stringpool.size(); j++)
        //     {
        //         if (i == j)
        //             continue;

        //         LinearExpr distance_j_x;

        //         for (int h = 1; h <= stringpool[i].size(); h++)
        //         {
        //             distance_j_x += CDEL * LinearExpr(yDel[StringIndelArc(i, h)]);

        //             for (int l = 1; l <= n; l++)
        //                 distance_j_x += CSUB * LinearExpr(yC[StringArc(i, h, l)]);
        //         }

        //         for (int l = 1; l <= n; l++)
        //             distance_j_x += CINS * LinearExpr(yIns[StringIndelArc(i, l)]);

        //         solver->MakeRowConstraint(dm(i, j) + distance_i_x >= distance_j_x);
        //     }

        //     for (int j = i + 1; j < stringpool.size(); j++)
        //     {
        //         LinearExpr distance_j_x;

        //         for (int h = 1; h <= stringpool[i].size(); h++)
        //         {
        //             distance_j_x += CDEL * LinearExpr(yDel[StringIndelArc(i, h)]);

        //             for (int l = 1; l <= n; l++)
        //                 distance_j_x += CSUB * LinearExpr(yC[StringArc(i, h, l)]);
        //         }

        //         for (int l = 1; l <= n; l++)
        //             distance_j_x += CINS * LinearExpr(yIns[StringIndelArc(i, l)]);

        //         solver->MakeRowConstraint(distance_i_x + distance_j_x >= dm(i, j));
        //     }
        // }

        // additional constraints
        LinearExpr sum_z;
        for (int k = 1; k <= n; k++)
            sum_z += LinearExpr(z[k]);

        solver->MakeRowConstraint(CINS * higher_str_size - LinearExpr(d) <=
                                  CINS * sum_z);

        solver->MakeRowConstraint(CDEL * sum_z <=
                                  CDEL * p + LinearExpr(d));

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
                LOG(INFO) << "The solver could not solve the problem. STATUS: " << result_status;
                return nullopt;
            }
        }

        vector<int> encoded_solution_str;
        string decoded_solution_str;

        encoded_solution_str.reserve(x.size());

        int i = 0;
        int solution_length = 0;
        for (MPVariable *var : z)
        {
            if (var->solution_value() < 0.5)
            {
                solution_length = i - 1;
                break;
            }

            i++;
        }

        i = 0;
        for (MPVariable *var : x)
        {
            if (i++ >= solution_length)
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
