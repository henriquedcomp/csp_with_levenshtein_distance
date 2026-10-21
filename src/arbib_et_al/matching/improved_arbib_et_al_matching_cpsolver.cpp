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
        optional<Solution> improved_arbib_et_al_matching(
            int n,
            const vector<vector<int>> &stringpool,
            const string_view &a,
            int p,
            const SymmetricMatrix<double> &dm,
            int higher_str_size,
            double solution_lower_bound)
        {
            CpModelBuilder cp_model;

            const double solution_ub = CSUB * p + CDEL * (higher_str_size - p);
            const Domain solution_domain(solution_lower_bound, solution_ub);
            const IntVar d = cp_model.NewIntVar(solution_domain);

            const Domain solution_string_domain(0, a.size() - 1);
            vector<IntVar> x(n);

            for (int i = 0; i < n; i++)
                x[i] = cp_model.NewIntVar(solution_string_domain);

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
                        for (int j = v; j <= n; j++)
                        {
                            LinearExpr sum = yF[StringArc(k, u, v)] +
                                             yC[StringArc(k, u, v)];

                            for (int h = 1; h < u; h++)
                            {
                                sum += yF[StringArc(k, h, j)] +
                                       yC[StringArc(k, h, j)];
                            }

                            cp_model.AddLessOrEqual(sum, 1.0);
                        }

                        for (int h = u; h <= stringpool[k].size(); h++)
                        {
                            LinearExpr sum = yF[StringArc(k, u, v)] +
                                             yC[StringArc(k, u, v)];

                            for (int j = 1; j < v; j++)
                            {
                                sum += yF[StringArc(k, h, j)] +
                                       yC[StringArc(k, h, j)];
                            }

                            cp_model.AddLessOrEqual(sum, 1.0);
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
                        cp_model.AddEquality(x[j - 1], lam).OnlyEnforceIf(yF[StringArc(k, h, j)]);
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

            // triangle inequality - it make results worse
            // for (int i = 0; i < stringpool.size(); i++)
            // {
            //     LinearExpr distance_i_x;

            //     for (int h = 1; h <= stringpool[i].size(); h++)
            //     {
            //         distance_i_x += CDEL * yDel[StringIndelArc(i, h)];

            //         for (int j = 1; j <= n; j++)
            //             distance_i_x += CSUB * yC[StringArc(i, h, j)];
            //     }

            //     for (int j = 1; j <= n; j++)
            //         distance_i_x += CINS * yIns[StringIndelArc(i, j)];

            //     for (int j = 0; j < stringpool.size(); j++)
            //     {
            //         if (i == j)
            //             continue;

            //         LinearExpr distance_j_x;

            //         for (int h = 1; h <= stringpool[i].size(); h++)
            //         {
            //             distance_j_x += CDEL * yDel[StringIndelArc(i, h)];

            //             for (int l = 1; l <= n; l++)
            //                 distance_j_x += CSUB * yC[StringArc(i, h, l)];
            //         }

            //         for (int l = 1; l <= n; l++)
            //             distance_j_x += CINS * yIns[StringIndelArc(i, l)];

            //         cp_model.AddGreaterOrEqual(dm(i, j) + distance_i_x, distance_j_x);
            //     }

            //     for (int j = i + 1; j < stringpool.size(); j++)
            //     {
            //         LinearExpr distance_j_x;

            //         for (int h = 1; h <= stringpool[i].size(); h++)
            //         {
            //             distance_j_x += CDEL * yDel[StringIndelArc(i, h)];

            //             for (int l = 1; l <= n; l++)
            //                 distance_j_x += CSUB * yC[StringArc(i, h, l)];
            //         }

            //         for (int l = 1; l <= n; l++)
            //             distance_j_x += CINS * yIns[StringIndelArc(i, l)];

            //         cp_model.AddGreaterOrEqual(distance_i_x + distance_j_x, dm(i, j));
            //     }
            // }

            // additional constraints
            LinearExpr sum_z;
            for (int k = 1; k <= n; k++)
                sum_z += z[k];

            cp_model.AddLessOrEqual(CINS * higher_str_size - d, CINS * sum_z);
            cp_model.AddLessOrEqual(CDEL * sum_z, CDEL * p + d);

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
            for (IntVar var : x)
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