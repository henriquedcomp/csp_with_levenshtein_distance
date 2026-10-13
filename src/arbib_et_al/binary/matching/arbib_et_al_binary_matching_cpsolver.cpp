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
        optional<Solution> arbib_et_al_binary_matching(
            int n,
            const vector<vector<int>> &stringpool,
            int p,
            int higher_str_size)
        {
            CpModelBuilder cp_model;

            const double solution_ub = CSUB * p + CDEL * (higher_str_size - p);
            const Domain solution_domain(0, solution_ub);
            const IntVar d = cp_model.NewIntVar(solution_domain);

            vector<BoolVar> x(n);

            for (int i = 0; i < n; i++)
                x[i] = cp_model.NewBoolVar();

            OffsetVector<BoolVar> z(1, n);
            for (int j = 1; j <= n; j++)
                z[j] = cp_model.NewBoolVar();

            for (int j = 1; j < n; j++)
                cp_model.AddLessOrEqual(z[j + 1], z[j]);

            cp_model.AddEquality(z[1], 1.0);

            map<StringArc, BoolVar> yF; // free substitutions
            map<StringArc, BoolVar> yC; // costly substitutions

            for (int k = 0; k < stringpool.size(); k++)
            {
                for (int h = 1; h <= stringpool[k].size(); h++)
                {
                    for (int j = 1; j <= n; j++)
                    {
                        yF[StringArc(k, h, j)] = cp_model.NewBoolVar();
                        yC[StringArc(k, h, j)] = cp_model.NewBoolVar();
                    }
                }
            }

            map<StringIndelArc, BoolVar> yDel;
            for (int k = 0; k < stringpool.size(); k++)
                for (int h = 1; h <= stringpool[k].size(); h++)
                    yDel[StringIndelArc(k, h)] = cp_model.NewBoolVar();

            map<StringIndelArc, BoolVar> yIns;
            for (int k = 0; k < stringpool.size(); k++)
                for (int j = 1; j <= n; j++)
                    yIns[StringIndelArc(k, j)] = cp_model.NewBoolVar();

            for (int k = 0; k < stringpool.size(); k++)
            {
                for (int h = 1; h <= stringpool[k].size(); h++)
                {
                    LinearExpr sum = yDel[StringIndelArc(k, h)];

                    for (int j = 1; j <= n; j++)
                        sum += yF[StringArc(k, h, j)] +
                               yC[StringArc(k, h, j)];

                    cp_model.AddEquality(sum, 1.0);
                }
            }

            for (int k = 0; k < stringpool.size(); k++)
            {
                for (int h = 1; h <= stringpool[k].size(); h++)
                {
                    LinearExpr sumJ, sumH, sumL;

                    for (int j = 1; j <= n; j++)
                        sumJ += yF[StringArc(k, h, j)] +
                                yC[StringArc(k, h, j)];

                    for (int u = h + 1; u <= stringpool[k].size(); u++)
                        sumH += yF[StringArc(k, u, 1)] +
                                yC[StringArc(k, u, 1)];

                    for (int u = 1; u < h; u++)
                        sumL += yF[StringArc(k, u, n)] +
                                yC[StringArc(k, u, n)];

                    cp_model.AddLessOrEqual(sumJ + sumH + sumL, 1.0);
                }
            }

            for (int k = 0; k < stringpool.size(); k++)
            {
                for (int j = 1; j <= n; j++)
                {
                    LinearExpr sum = yIns[StringIndelArc(k, j)];

                    for (int h = 1; h <= stringpool[k].size(); h++)
                        sum += yF[StringArc(k, h, j)] +
                               yC[StringArc(k, h, j)];

                    cp_model.AddEquality(sum, z[j]);
                }
            }

            for (int k = 0; k < stringpool.size(); k++)
            {
                for (int j = 1; j <= n; j++)
                {
                    LinearExpr sum_h, sumH, sumL;

                    for (int h = 1; h <= stringpool[k].size(); h++)
                        sum_h += yF[StringArc(k, h, j)] +
                                 yC[StringArc(k, h, j)];

                    for (int v = j + 1; v <= n; v++)
                        sumH += yF[StringArc(k, 1, v)] +
                                yC[StringArc(k, 1, v)];

                    for (int v = 1; v < j; v++)
                        sumL += yF[StringArc(k, stringpool[k].size(), v)] +
                                yC[StringArc(k, stringpool[k].size(), v)];

                    cp_model.AddLessOrEqual(sum_h + sumH + sumL, 1.0);
                }
            }

            // nocross constraint
            for (int k = 0; k < stringpool.size(); k++)
            {
                for (int u = 1; u <= stringpool[k].size(); u++)
                {
                    for (int v = 1; v <= n; v++)
                    {
                        LinearExpr lhs;

                        for (int j = 1; j <= n; j++)
                        {
                            if (j == v)
                                continue;

                            lhs += yF[StringArc(k, u, j)] +
                                   yC[StringArc(k, u, j)];
                        }

                        for (int h = 1; h <= stringpool[k].size(); h++)
                        {
                            if (h == u)
                                continue;

                            lhs += yF[StringArc(k, h, v)] +
                                   yC[StringArc(k, h, v)];
                        }

                        for (int h = 1; h < u; h++)
                        {
                            for (int j = v + 1; j <= n; j++)
                            {
                                lhs += yF[StringArc(k, h, j)] +
                                       yC[StringArc(k, h, j)];
                            }
                        }

                        for (int h = u + 1; h <= stringpool[k].size(); h++)
                        {
                            for (int j = 1; j < v; j++)
                            {
                                lhs += yF[StringArc(k, h, j)] +
                                       yC[StringArc(k, h, j)];
                            }
                        }

                        cp_model.AddLessOrEqual(lhs + yDel[StringIndelArc(k, u)],
                                                stringpool[k].size() * (1.0 - yF[StringArc(k, u, v)]) - yC[StringArc(k, u, v)]);
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
                        if (lam == 1)
                            cp_model.AddLessOrEqual(yF[StringArc(k, h, j)] - x[j - 1], 0.0);
                        else
                            cp_model.AddLessOrEqual(yF[StringArc(k, h, j)] + x[j - 1], 1.0);
                    }
                }
            }

            for (int k = 0; k < stringpool.size(); k++)
            {
                LinearExpr distance;

                for (int h = 1; h <= stringpool[k].size(); h++)
                {
                    distance += CDEL * yDel[StringIndelArc(k, h)];

                    for (int j = 1; j <= n; j++)
                        distance += CSUB * yC[StringArc(k, h, j)];
                }

                for (int j = 1; j <= n; j++)
                    distance += CINS * yIns[StringIndelArc(k, j)];

                cp_model.AddLessOrEqual(distance, d);
            }

            cp_model.Minimize(d);

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

            encoded_solution_str.reserve(x.size());

            int i = 0;
            int solution_length = 0;
            for (BoolVar var : z)
            {
                if (SolutionIntegerValue(response, var) < 0.5)
                {
                    solution_length = i - 1;
                    break;
                }

                i++;
            }

            i = 0;
            for (BoolVar var : x)
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
                static_cast<double>(SolutionIntegerValue(response, d)),
                status,
                response.wall_time()}; // seconds

            return solution;
        }
    }
}