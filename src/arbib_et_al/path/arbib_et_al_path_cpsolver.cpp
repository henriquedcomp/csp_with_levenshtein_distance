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
        optional<Solution> arbib_et_al_path(
            int sol_size_ub,
            int sol_size_lb,
            const vector<vector<int>> &input_strings,
            const string_view &alphabet,
            double sol_ub,
            double sol_lb)
        {
            CpModelBuilder cp_model;

            const Domain solution_domain(sol_lb, sol_ub);
            const IntVar r = cp_model.NewIntVar(solution_domain);

            const Domain solution_string_domain(0, alphabet.size() - 1);
            vector<IntVar> t(sol_size_ub);

            for (int i = 0; i < sol_size_ub; i++)
                t[i] = cp_model.NewIntVar(solution_string_domain);

            OffsetVector<BoolVar> l(sol_size_lb, sol_size_ub + 1);
            for (int k = sol_size_lb; k <= sol_size_ub + 1; k++)
                l[k] = cp_model.NewBoolVar();

            for (int k = sol_size_lb + 1; k < sol_size_ub; k++)
                cp_model.AddGreaterOrEqual(l[k], l[k + 1]); // (4.35)

            cp_model.AddEquality(l[sol_size_lb], 1.0); // (4.36)
            cp_model.AddEquality(l[sol_size_ub + 1], 0.0); // (4.37)

            for (int i = 0; i < input_strings.size(); i++)
            {
                int nk = input_strings[i].size();

                vector<Arc> arcs;
                GraphConstruction::generate_arcs(nk, sol_size_ub, arcs);

                vector<Vertex> nodes;
                GraphConstruction::generate_nodes(arcs, nodes);

                vector<FreeSlantArc> free_slant_arcs;
                GraphConstruction::generate_free_slant_arcs(input_strings[i], nk, sol_size_ub, free_slant_arcs);

                map<Arc, BoolVar> flow;
                for (const auto &arc : arcs)
                    flow[arc] = cp_model.NewBoolVar();

                map<FreeSlantArc, BoolVar> zero_cost_flow;
                for (const auto &arc : free_slant_arcs)
                    zero_cost_flow[arc] = cp_model.NewBoolVar();

                for (const auto &arc : arcs)
                {
                    int j = get<1>(get<1>(arc));

                    if (sol_size_lb < j && j <= sol_size_ub)
                        cp_model.AddLessOrEqual(flow[arc], l[j]); // (4.39)
                }

                for (const auto &arc : free_slant_arcs)
                {
                    int j = get<1>(get<1>(arc));

                    if (sol_size_lb < j && j <= sol_size_ub)
                        cp_model.AddLessOrEqual(zero_cost_flow[arc], l[j]); // (4.39)
                }

                for (const auto &zcf : zero_cost_flow)
                    cp_model.AddEquality(t[get<1>(get<1>(zcf.first)) - 1], get<2>(zcf.first)).OnlyEnforceIf(zcf.second); // (4.49)

                LinearExpr source;

                for (const auto &arc : arcs)
                    if (get<0>(arc) == make_tuple(0, 0))
                        source += flow[arc];

                for (const auto &arc : free_slant_arcs)
                    if (get<0>(arc) == make_tuple(0, 0))
                        source += zero_cost_flow[arc];

                cp_model.AddEquality(source, 1.0); // (4.31)

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
                            outflow += flow[arc];

                        if (get<1>(arc) == vertex)
                            inflow += flow[arc];
                    }

                    for (const auto &arc : free_slant_arcs)
                    {
                        if (get<0>(arc) == vertex)
                            outflow += zero_cost_flow[arc];

                        if (get<1>(arc) == vertex)
                            inflow += zero_cost_flow[arc];
                    }

                    cp_model.AddEquality(outflow, inflow); // (4.32)
                }

                for (int k = sol_size_lb; k <= sol_size_ub; k++)
                {
                    LinearExpr outflow, inflow;
                    Vertex node(nk, k);

                    for (const auto &arc : arcs)
                    {
                        if (get<0>(arc) == node)
                            outflow += flow[arc];

                        if (get<1>(arc) == node)
                            inflow += flow[arc];
                    }

                    for (const auto &arc : free_slant_arcs)
                    {
                        if (get<0>(arc) == node)
                            outflow += zero_cost_flow[arc];

                        if (get<1>(arc) == node)
                            inflow += zero_cost_flow[arc];
                    }

                    cp_model.AddEquality(inflow - outflow, l[k] - l[k + 1]); // (4.38)
                }

                LinearExpr distance;

                for (const auto &value : flow)
                    distance += get<2>(value.first) * value.second;

                cp_model.AddLessOrEqual(distance, r); // (4.30)
            }

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
                    solution_length = sol_size_lb + i - 1;
                    break;
                }

                i++;
            }

            i = 0;
            for (IntVar var : t)
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