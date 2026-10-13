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
            int m,
            const vector<vector<int>> &s,
            const string_view &a,
            const SymmetricMatrix<double> &dm,
            int lower_str_size,
            int higher_str_size,
            double solution_lower_bound)
        {
            CpModelBuilder cp_model;

            const double solution_ub = CSUB * lower_str_size + CDEL * (higher_str_size - lower_str_size);
            const Domain solution_domain(solution_lower_bound, solution_ub);
            const Domain solution_size_domain(0, m);

            const IntVar d = cp_model.NewIntVar(solution_domain);
            const IntVar l = cp_model.NewIntVar(solution_size_domain);

            const Domain alphabet_domain(0, a.size() - 1);
            vector<IntVar> t(m);

            for (int i = 0; i < m; i++)
                t[i] = cp_model.NewIntVar(alphabet_domain);

            vector<vector<vector<BoolVar>>> x(s.size() + 1);
            vector<vector<vector<BoolVar>>> y(s.size() + 1);
            vector<vector<vector<BoolVar>>> z(s.size() + 1);
            vector<vector<vector<BoolVar>>> g(s.size() + 1);
            vector<vector<vector<BoolVar>>> h(s.size() + 1);

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

            LinearExpr diff = m - l;

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

                cp_model.AddLessOrEqual(distance, d);
                cp_model.AddExactlyOne({x[k][1][0], y[k][0][1], z[k][1][1]});
                cp_model.AddEquality(x[k][nk][0], y[k][nk][1]);
                cp_model.AddEquality(x[k][1][m], y[k][0][m]);
                cp_model.AddExactlyOne({x[k][nk][m], y[k][nk][m], z[k][nk][m]});

                for (int i = 1; i < nk; i++)
                {
                    cp_model.AddEquality(x[k][i][0], x[k][i + 1][0] + y[k][i][1] + z[k][i + 1][1]);

                    cp_model.AddEquality(x[k][i][m] + y[k][i][m] + z[k][i][m], x[k][i + 1][m]);

                    for (int j = 1; j < m; j++)
                    {
                        cp_model.AddEquality(
                            x[k][i][j] + y[k][i][j] + z[k][i][j],
                            x[k][i + 1][j] + y[k][i][j + 1] + z[k][i + 1][j + 1]);
                    }
                }

                for (int j = 1; j < m; j++)
                {
                    cp_model.AddEquality(y[k][0][j], x[k][1][j] + y[k][0][j + 1] + z[k][1][j + 1]);
                    cp_model.AddEquality(x[k][nk][j] + y[k][nk][j] + z[k][nk][j], y[k][nk][j + 1]);
                }

                for (int j = 1; j <= m; j++)
                {
                    cp_model.AddLessOrEqual(j, l).OnlyEnforceIf(Not(y[k][nk][j]));

                    for (int i = 1; i <= nk; i++)
                    {
                        cp_model.AddEquality(s[k][i - 1], t[j - 1]).OnlyEnforceIf(Not(g[k][i][j]));
                        cp_model.AddEquality(h[k][i][j], 1).OnlyEnforceIf({z[k][i][j], g[k][i][j]});
                    }
                }
            }

            // triangle inequality - it makes result worse
            // for (int k1 = 0; k1 < s.size(); k1++)
            // {
            //     LinearExpr distance_k1_t;
            //     size_t nk1 = s[k1].size();

            //     for (int i = 1; i <= nk1; i++)
            //     {
            //         distance_k1_t += CDEL * x[k1][i][0];

            //         for (int j = 1; j <= m; j++)
            //         {
            //             distance_k1_t += CDEL * x[k1][i][j] +
            //                              CINS * y[k1][i][j] +
            //                              CSUB * h[k1][i][j];
            //         }
            //     }

            //     for (int j = 1; j <= m; j++)
            //         distance_k1_t += CINS * y[k1][0][j];

            //     distance_k1_t -= CINS * diff;

            //     for (int k2 = 0; k2 < s.size(); k2++)
            //     {
            //         if (k2 == k1)
            //             continue;

            //         LinearExpr distance_k2_t;

            //         size_t nk2 = s[k2].size();

            //         for (int i = 1; i <= nk2; i++)
            //         {
            //             distance_k2_t += CDEL * x[k2][i][0];

            //             for (int j = 1; j <= m; j++)
            //             {
            //                 distance_k2_t += CDEL * x[k2][i][j] +
            //                                  CINS * y[k2][i][j] +
            //                                  CSUB * h[k2][i][j];
            //             }
            //         }

            //         for (int j = 1; j <= m; j++)
            //             distance_k2_t += CINS * y[k2][0][j];

            //         distance_k2_t -= CINS * diff;

            //         cp_model.AddGreaterOrEqual(dm(k1, k2) + distance_k1_t, distance_k2_t);
            //     }

            //     for (int k2 = k1 + 1; k2 < s.size(); k2++)
            //     {
            //         LinearExpr distance_k2_t;

            //         size_t nk2 = s[k2].size();

            //         for (int i = 1; i <= nk2; i++)
            //         {
            //             distance_k2_t += CDEL * x[k2][i][0];

            //             for (int j = 1; j <= m; j++)
            //             {
            //                 distance_k2_t += CDEL * x[k2][i][j] +
            //                                  CINS * y[k2][i][j] +
            //                                  CSUB * h[k2][i][j];
            //             }
            //         }

            //         for (int j = 1; j <= m; j++)
            //             distance_k2_t += CINS * y[k2][0][j];

            //         distance_k2_t -= CINS * diff;

            //         cp_model.AddGreaterOrEqual(distance_k1_t + distance_k2_t, dm(k1, k2));
            //     }
            // }

            // additional constraints
            cp_model.AddLessOrEqual(CINS * higher_str_size - d, l * CINS);
            cp_model.AddLessOrEqual(l * CDEL, lower_str_size * CDEL + d);

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
                static_cast<double>(SolutionIntegerValue(response, d)),
                status,
                response.wall_time()}; // seconds

            return solution;
        }
    }
}