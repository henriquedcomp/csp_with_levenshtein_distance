#include "formulations.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "offset_vector.hpp"
#include "graph_construction.hpp"
#include "ortools/sat/cp_model.h"

using namespace std;

using StringArc = tuple<int, int, int>;
using StringIndelArc = tuple<int, int>;

namespace operations_research
{
    namespace sat
    {
        optional<Solution> improved_arbib_et_al_binary_matching(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            int higher_str_size,
            int lower_str_size,
            double sol_ub,
            double sol_lb)
        {
            CpModelBuilder cp_model;

            const Domain solution_domain(sol_lb, sol_ub);
            const IntVar r = cp_model.NewIntVar(solution_domain);

            vector<BoolVar> t(sol_size_ub);

            for (int i = 0; i < sol_size_ub; i++)
                t[i] = cp_model.NewBoolVar();

            OffsetVector<BoolVar> l(sol_size_lb + 1, sol_size_ub);
            for (int j = sol_size_lb + 1; j <= sol_size_ub; j++)
                l[j] = cp_model.NewBoolVar();

            for (int j = sol_size_lb + 1; j < sol_size_ub; j++)
                cp_model.AddLessOrEqual(l[j + 1], l[j]); // (4.58)

            map<StringArc, BoolVar> yF; // free substitutions
            map<StringArc, BoolVar> yC; // costly substitutions

            for (int k = 0; k < input_strings.size(); k++)
            {
                for (int h = 1; h <= input_strings[k].size(); h++)
                {
                    for (int j = 1; j <= sol_size_ub; j++)
                    {
                        yF[StringArc(k, h, j)] = cp_model.NewBoolVar();
                        yC[StringArc(k, h, j)] = cp_model.NewBoolVar();
                    }
                }
            }

            map<StringIndelArc, BoolVar> yDel;
            for (int k = 0; k < input_strings.size(); k++)
                for (int h = 1; h <= input_strings[k].size(); h++)
                    yDel[StringIndelArc(k, h)] = cp_model.NewBoolVar();

            map<StringIndelArc, BoolVar> yIns;
            for (int k = 0; k < input_strings.size(); k++)
                for (int j = 1; j <= sol_size_ub; j++)
                    yIns[StringIndelArc(k, j)] = cp_model.NewBoolVar();

            for (int k = 0; k < input_strings.size(); k++)
            {
                for (int h = 1; h <= input_strings[k].size(); h++)
                {
                    LinearExpr sum = yDel[StringIndelArc(k, h)];

                    for (int j = 1; j <= sol_size_ub; j++)
                        sum += yF[StringArc(k, h, j)] +
                               yC[StringArc(k, h, j)];

                    cp_model.AddEquality(sum, 1.0); // (4.52)
                }
            }

            for (int k = 0; k < input_strings.size(); k++)
            {
                for (int h = 1; h <= input_strings[k].size(); h++)
                {
                    LinearExpr sumJ, sumH, sumL;

                    for (int j = 1; j <= sol_size_ub; j++)
                        sumJ += yF[StringArc(k, h, j)] +
                                yC[StringArc(k, h, j)];

                    for (int u = h + 1; u <= input_strings[k].size(); u++)
                        sumH += yF[StringArc(k, u, 1)] +
                                yC[StringArc(k, u, 1)];

                    for (int u = 1; u < h; u++)
                        sumL += yF[StringArc(k, u, sol_size_ub)] +
                                yC[StringArc(k, u, sol_size_ub)];

                    cp_model.AddLessOrEqual(sumJ + sumH + sumL, 1.0); // (4.59)
                }
            }

            for (int k = 0; k < input_strings.size(); k++)
            {
                for (int j = 1; j <= sol_size_lb; j++)
                {
                    LinearExpr sum = yIns[StringIndelArc(k, j)];

                    for (int h = 1; h <= input_strings[k].size(); h++)
                        sum += yF[StringArc(k, h, j)] +
                               yC[StringArc(k, h, j)];

                    cp_model.AddEquality(sum, 1.0); // (4.53)
                }

                for (int j = sol_size_lb + 1; j <= sol_size_ub; j++)
                {
                    LinearExpr sum = yIns[StringIndelArc(k, j)];

                    for (int h = 1; h <= input_strings[k].size(); h++)
                        sum += yF[StringArc(k, h, j)] +
                               yC[StringArc(k, h, j)];

                    cp_model.AddEquality(sum, l[j]); // (4.54)
                }
            }

            for (int k = 0; k < input_strings.size(); k++)
            {
                for (int j = 1; j <= sol_size_ub; j++)
                {
                    LinearExpr sum_h, sumH, sumL;

                    for (int h = 1; h <= input_strings[k].size(); h++)
                        sum_h += yF[StringArc(k, h, j)] +
                                 yC[StringArc(k, h, j)];

                    for (int v = j + 1; v <= sol_size_ub; v++)
                        sumH += yF[StringArc(k, 1, v)] +
                                yC[StringArc(k, 1, v)];

                    for (int v = 1; v < j; v++)
                        sumL += yF[StringArc(k, input_strings[k].size(), v)] +
                                yC[StringArc(k, input_strings[k].size(), v)];

                    cp_model.AddLessOrEqual(sum_h + sumH + sumL, 1.0); // (4.60)
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
                            LinearExpr sum = yF[StringArc(k, u, v)] +
                                             yC[StringArc(k, u, v)];

                            for (int h = 1; h < u; h++)
                            {
                                sum += yF[StringArc(k, h, j)] +
                                       yC[StringArc(k, h, j)];
                            }

                            cp_model.AddLessOrEqual(sum, 1.0); // (4.64)
                        }

                        for (int h = u + 1; h <= input_strings[k].size(); h++)
                        {
                            LinearExpr sum = yF[StringArc(k, u, v)] +
                                             yC[StringArc(k, u, v)];

                            for (int j = 1; j < v; j++)
                            {
                                sum += yF[StringArc(k, h, j)] +
                                       yC[StringArc(k, h, j)];
                            }

                            cp_model.AddLessOrEqual(sum, 1.0); // (4.65)
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
                        if (lam == 1)
                            cp_model.AddLessOrEqual(yF[StringArc(k, h, j)] - t[j - 1], 0.0); // (4.56)
                        else
                            cp_model.AddLessOrEqual(yF[StringArc(k, h, j)] + t[j - 1], 1.0); // (4.57)
                    }
                }
            }

            for (int k = 0; k < input_strings.size(); k++)
            {
                LinearExpr distance;

                for (int h = 1; h <= input_strings[k].size(); h++)
                {
                    distance += CDEL * yDel[StringIndelArc(k, h)];

                    for (int j = 1; j <= sol_size_ub; j++)
                        distance += CSUB * yC[StringArc(k, h, j)];
                }

                for (int j = 1; j <= sol_size_ub; j++)
                    distance += CINS * yIns[StringIndelArc(k, j)];

                cp_model.AddLessOrEqual(distance, r); // (4.51)
            }

            // additional constraints
            LinearExpr sum_z = sol_size_lb;
            for (int k = sol_size_lb + 1; k <= sol_size_ub; k++)
                sum_z += l[k];

            cp_model.AddLessOrEqual(CINS * higher_str_size - r, CINS * sum_z); // (4.43)
            cp_model.AddLessOrEqual(CDEL * sum_z, CDEL * lower_str_size + r);  // (4.44)

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
            int solution_length = sol_size_lb;
            for (BoolVar var : l)
            {
                if (SolutionIntegerValue(response, var) < 0.5)
                {
                    solution_length = sol_size_lb + i;
                    break;
                }

                i++;
            }

            i = 0;
            for (BoolVar var : t)
            {
                if (i++ >= solution_length)
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
                response.wall_time()}; // seconds

            return solution;
        }
    }
}