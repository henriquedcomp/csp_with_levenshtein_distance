#include "formulations.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "offset_vector.hpp"
#include "graph_construction.hpp"

#include "ortools/sat/cp_model.h"

using namespace std;

namespace operations_research
{
    namespace sat
    {
        optional<Solution> improved_arbib_et_al_binary_path(
            int n,
            const vector<vector<int>> &stringpool,
            int p,
            const SymmetricMatrix<double> &dm,
            int higher_str_size,
            double solution_lower_bound)
        {
            CpModelBuilder cp_model;

            const double solution_ub = CSUB * p + CDEL * (higher_str_size - p);
            const Domain solution_domain(solution_lower_bound, solution_ub);
            const IntVar d = cp_model.NewIntVar(solution_domain);

            vector<BoolVar> x(n);

            for (int i = 0; i < n; i++)
                x[i] = cp_model.NewBoolVar();

            OffsetVector<BoolVar> u(0, n + 1);
            for (int k = 0; k <= n + 1; k++)
                u[k] = cp_model.NewBoolVar();

            for (int k = 1; k < n; k++)
                cp_model.AddGreaterOrEqual(u[k], u[k + 1]);

            cp_model.AddEquality(u[0], 1.0);
            cp_model.AddEquality(u[n + 1], 0.0);

            vector<map<Arc, BoolVar>> flow(stringpool.size());

            for (int i = 0; i < stringpool.size(); i++)
            {
                int nk = stringpool[i].size();

                vector<Arc> arcs;
                GraphConstruction::generate_arcs(nk, n, arcs);

                vector<Vertex> nodes;
                GraphConstruction::generate_nodes(arcs, nodes);

                vector<FreeSlantArc> free_slant_arcs;
                GraphConstruction::generate_free_slant_arcs(stringpool[i], nk, n, free_slant_arcs);

                auto &internal_flow = flow[i];
                for (const auto &arc : arcs)
                    internal_flow[arc] = cp_model.NewBoolVar();

                map<FreeSlantArc, BoolVar> zero_cost_flow;
                for (const auto &arc : free_slant_arcs)
                    zero_cost_flow[arc] = cp_model.NewBoolVar();

                for (const auto &arc : arcs)
                {
                    int j = get<1>(get<1>(arc));

                    if (0 < j && j <= n)
                        cp_model.AddLessOrEqual(internal_flow[arc], u[j]);
                }

                for (const auto &arc : free_slant_arcs)
                {
                    int j = get<1>(get<1>(arc));

                    if (0 < j && j <= n)
                        cp_model.AddLessOrEqual(zero_cost_flow[arc], u[j]);
                }

                for (const auto &zcf : zero_cost_flow)
                {
                    if (get<2>(zcf.first) == 0)
                    {
                        cp_model.AddLessOrEqual(zcf.second, 1 - x[get<1>(get<1>(zcf.first)) - 1]);
                    }
                    else
                    {
                        cp_model.AddLessOrEqual(zcf.second, x[get<1>(get<1>(zcf.first)) - 1]);
                    }
                }

                LinearExpr source;

                for (const auto &arc : arcs)
                    if (get<0>(arc) == make_tuple(0, 0))
                        source += internal_flow[arc];

                for (const auto &arc : free_slant_arcs)
                    if (get<0>(arc) == make_tuple(0, 0))
                        source += zero_cost_flow[arc];

                cp_model.AddEquality(source, 1.0);

                for (const auto &vertex : nodes)
                {
                    int i = get<0>(vertex);
                    int j = get<1>(vertex);

                    if (vertex == make_tuple(0, 0) || i == nk)
                        continue;

                    LinearExpr outflow, inflow;

                    for (const auto &arc : arcs)
                    {
                        if (get<0>(arc) == vertex)
                            outflow += internal_flow[arc];

                        if (get<1>(arc) == vertex)
                            inflow += internal_flow[arc];
                    }

                    for (const auto &arc : free_slant_arcs)
                    {
                        if (get<0>(arc) == vertex)
                            outflow += zero_cost_flow[arc];

                        if (get<1>(arc) == vertex)
                            inflow += zero_cost_flow[arc];
                    }

                    cp_model.AddEquality(outflow, inflow);
                }

                for (int k = 0; k <= n; k++)
                {
                    LinearExpr outflow, inflow;
                    Vertex node(nk, k);

                    for (const auto &arc : arcs)
                    {
                        if (get<0>(arc) == node)
                            outflow += internal_flow[arc];

                        if (get<1>(arc) == node)
                            inflow += internal_flow[arc];
                    }

                    for (const auto &arc : free_slant_arcs)
                    {
                        if (get<0>(arc) == node)
                            outflow += zero_cost_flow[arc];

                        if (get<1>(arc) == node)
                            inflow += zero_cost_flow[arc];
                    }

                    cp_model.AddEquality(inflow - outflow, u[k] - u[k + 1]);
                }

                LinearExpr distance;

                for (const auto &value : internal_flow)
                    distance += get<2>(value.first) * value.second;

                cp_model.AddLessOrEqual(distance, d);
            }

            // triangle inequality - it makes result worse
            // for (int i = 0; i < stringpool.size(); i++)
            // {
            //     LinearExpr distance_i_x;

            //     for (const auto &value : flow[i])
            //         distance_i_x += value.second;

            //     for (int j = 0; j < stringpool.size(); j++)
            //     {
            //         if (i == j)
            //             continue;

            //         LinearExpr distance_j_x;

            //         for (const auto &value : flow[j])
            //             distance_j_x += value.second;

            //         cp_model.AddGreaterOrEqual(dm(i, j) + distance_i_x, distance_j_x);
            //     }

            //     for (int j = i + 1; j < stringpool.size(); j++)
            //     {
            //         LinearExpr distance_j_x;

            //         for (const auto &value : flow[j])
            //             distance_j_x += value.second;

            //         cp_model.AddGreaterOrEqual(distance_i_x + distance_j_x, dm(i, j));
            //     }
            // }

            // additional constraints
            LinearExpr sum_u;
            for (int k = 1; k <= n; k++)
                sum_u += u[k];

            cp_model.AddLessOrEqual(CINS * higher_str_size - d, CINS * sum_u);
            cp_model.AddLessOrEqual(CDEL * sum_u, CDEL * p + d);

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
            for (BoolVar var : u)
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