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
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
        int higher_str_size,
        int lower_str_size,
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

        vector<MPVariable *> t;
        solver->MakeNumVarArray(sol_size_ub, 0, alphabet.size() - 1, "t", &t);

        OffsetVector<MPVariable *> l(sol_size_lb + 1, sol_size_ub);
        for (int j = sol_size_lb + 1; j <= sol_size_ub; j++)
        {
            string variable_name = format("l_{}", j);
            l[j] = solver->MakeBoolVar(variable_name);
        }

        for (int j = sol_size_lb + 1; j < sol_size_ub; j++)
            solver->MakeRowConstraint(LinearExpr(l[j + 1]) <= LinearExpr(l[j])); // (4.58)

        map<StringArc, MPVariable *> yF; // free substitutions
        map<StringArc, MPVariable *> yC; // costly substitutions

        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int h = 1; h <= input_strings[k].size(); h++)
            {
                for (int j = 1; j <= sol_size_ub; j++)
                {
                    string free_variable_name = format("yF_{}_{}_{}", k, h, j);
                    string costly_variable_name = format("yC_{}_{}_{}", k, h, j);

                    yF[StringArc(k, h, j)] = solver->MakeBoolVar(free_variable_name);
                    yC[StringArc(k, h, j)] = solver->MakeBoolVar(costly_variable_name);
                }
            }
        }

        map<StringIndelArc, MPVariable *> yDel;
        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int h = 1; h <= input_strings[k].size(); h++)
            {
                string variable_name = format("yDel_{}_{}", k, h);

                yDel[StringIndelArc(k, h)] = solver->MakeBoolVar(variable_name);
            }
        }

        map<StringIndelArc, MPVariable *> yIns;
        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int j = 1; j <= sol_size_ub; j++)
            {
                string variable_name = format("yIns_{}_{}", k, j);

                yIns[StringIndelArc(k, j)] = solver->MakeBoolVar(variable_name);
            }
        }

        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int h = 1; h <= input_strings[k].size(); h++)
            {
                LinearExpr sum = LinearExpr(yDel[StringIndelArc(k, h)]);

                for (int j = 1; j <= sol_size_ub; j++)
                    sum += LinearExpr(yF[StringArc(k, h, j)]) +
                           LinearExpr(yC[StringArc(k, h, j)]);

                solver->MakeRowConstraint(sum == 1.0); // (4.52)
            }
        }

        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int h = 1; h <= input_strings[k].size(); h++)
            {
                LinearExpr sumJ, sumH, sumL;

                for (int j = 1; j <= sol_size_ub; j++)
                    sumJ += LinearExpr(yF[StringArc(k, h, j)]) +
                            LinearExpr(yC[StringArc(k, h, j)]);

                for (int u = h + 1; u <= input_strings[k].size(); u++)
                    sumH += LinearExpr(yF[StringArc(k, u, 1)]) +
                            LinearExpr(yC[StringArc(k, u, 1)]);

                for (int u = 1; u < h; u++)
                    sumL += LinearExpr(yF[StringArc(k, u, sol_size_ub)]) +
                            LinearExpr(yC[StringArc(k, u, sol_size_ub)]);

                solver->MakeRowConstraint(sumJ + sumH + sumL <= 1.0); // (4.59)
            }
        }

        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int j = 1; j <= sol_size_lb; j++)
            {
                LinearExpr sum = LinearExpr(yIns[StringIndelArc(k, j)]);

                for (int h = 1; h <= input_strings[k].size(); h++)
                    sum += LinearExpr(yF[StringArc(k, h, j)]) +
                           LinearExpr(yC[StringArc(k, h, j)]);

                solver->MakeRowConstraint(sum == 1.0); // (4.53)
            }

            for (int j = sol_size_lb + 1; j <= sol_size_ub; j++)
            {
                LinearExpr sum = LinearExpr(yIns[StringIndelArc(k, j)]);

                for (int h = 1; h <= input_strings[k].size(); h++)
                    sum += LinearExpr(yF[StringArc(k, h, j)]) +
                           LinearExpr(yC[StringArc(k, h, j)]);

                solver->MakeRowConstraint(sum == l[j]); // (4.54)
            }
        }

        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int j = 1; j <= sol_size_ub; j++)
            {
                LinearExpr sum_h, sumH, sumL;

                for (int h = 1; h <= input_strings[k].size(); h++)
                    sum_h += LinearExpr(yF[StringArc(k, h, j)]) +
                             LinearExpr(yC[StringArc(k, h, j)]);

                for (int v = j + 1; v <= sol_size_ub; v++)
                    sumH += LinearExpr(yF[StringArc(k, 1, v)]) +
                            LinearExpr(yC[StringArc(k, 1, v)]);

                for (int v = 1; v < j; v++)
                    sumL += LinearExpr(yF[StringArc(k, input_strings[k].size(), v)]) +
                            LinearExpr(yC[StringArc(k, input_strings[k].size(), v)]);

                solver->MakeRowConstraint(sum_h + sumH + sumL <= 1.0); // (4.60)
            }
        }

        // nocross constraint reimagined
        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int u = 1; u <= input_strings[k].size(); u++)
            {
                for (int v = 1; v <= sol_size_ub; v++)
                {
                    for (int j = v + 1; j <= sol_size_ub; j++)
                    {
                        LinearExpr sum = LinearExpr(yF[StringArc(k, u, v)]) +
                                         LinearExpr(yC[StringArc(k, u, v)]);

                        for (int h = 1; h < u; h++)
                        {
                            sum += LinearExpr(yF[StringArc(k, h, j)]) +
                                   LinearExpr(yC[StringArc(k, h, j)]);
                        }

                        solver->MakeRowConstraint(sum <= 1.0); // (4.64)
                    }

                    for (int h = u + 1; h <= input_strings[k].size(); h++)
                    {
                        LinearExpr sum = LinearExpr(yF[StringArc(k, u, v)]) +
                                         LinearExpr(yC[StringArc(k, u, v)]);

                        for (int j = 1; j < v; j++)
                        {
                            sum += LinearExpr(yF[StringArc(k, h, j)]) +
                                   LinearExpr(yC[StringArc(k, h, j)]);
                        }

                        solver->MakeRowConstraint(sum <= 1.0); // (4.65)
                    }
                }
            }
        }

        for (int k = 0; k < input_strings.size(); k++)
        {
            for (int h = 1; h <= input_strings[k].size(); h++)
            {
                int lam = input_strings[k][h - 1];

                for (int j = 1; j <= sol_size_ub; j++)
                {
                    string name = format("e_{}_{}_{}", k, h, j);
                    MPVariable *e = solver->MakeBoolVar(name);

                    solver->MakeRowConstraint(lam - LinearExpr(t[j - 1]) <= alphabet.size() * LinearExpr(e)); // (4.43)

                    solver->MakeRowConstraint(LinearExpr(t[j - 1]) - lam <= alphabet.size() * LinearExpr(e)); // (4.44)

                    solver->MakeRowConstraint(LinearExpr(yF[StringArc(k, h, j)]) + LinearExpr(e) <= 1.0); // (4.45)
                }
            }
        }

        for (int k = 0; k < input_strings.size(); k++)
        {
            LinearExpr distance;

            for (int h = 1; h <= input_strings[k].size(); h++)
            {
                distance += CDEL * LinearExpr(yDel[StringIndelArc(k, h)]);

                for (int j = 1; j <= sol_size_ub; j++)
                    distance += CSUB * LinearExpr(yC[StringArc(k, h, j)]);
            }

            for (int j = 1; j <= sol_size_ub; j++)
                distance += CINS * LinearExpr(yIns[StringIndelArc(k, j)]);

            solver->MakeRowConstraint(distance <= r); // (4.51)
        }

        // additional constraints
        LinearExpr sum_z = sol_size_lb;
        for (int k = sol_size_lb + 1; k <= sol_size_ub; k++)
            sum_z += LinearExpr(l[k]);

        solver->MakeRowConstraint(CINS * higher_str_size - LinearExpr(r) <=
                                  CINS * sum_z); // (4.43)

        solver->MakeRowConstraint(CDEL * sum_z <=
                                  CDEL * lower_str_size + LinearExpr(r)); // (4.44)

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
                LOG(INFO) << "The solver could not solve the problem. STATUS: " << result_status;
                return nullopt;
            }
        }

        vector<int> encoded_solution_str;
        string decoded_solution_str;

        encoded_solution_str.reserve(t.size());

        int i = 0;
        int solution_length = sol_size_lb;
        for (MPVariable *var : l)
        {
            if (var->solution_value() < 0.5)
            {
                solution_length = sol_size_lb + i;
                break;
            }

            i++;
        }

        i = 0;
        for (MPVariable *var : t)
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
            r->solution_value(),
            status,
            absl::ToDoubleSeconds(solver->DurationSinceConstruction())};

        return solution;
    }
}
