#include <memory>
#include <map>

#include "formulations.hpp"
#include "string_codec.hpp"
#include "config.hpp"
#include "offset_vector.hpp"
#include "graph_construction.hpp"

#include "ortools/linear_solver/linear_solver.h"

using namespace std;

namespace operations_research
{
    optional<Solution> improved_arbib_et_al_binary_path(
        int n,
        const vector<vector<int>> &stringpool,
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
        solver->MakeBoolVarArray(n, "x", &x);

        OffsetVector<MPVariable *> u(0, n + 1);
        for (int k = 0; k <= n + 1; k++)
        {
            string variable_name = format("u_{}", k);
            u[k] = solver->MakeBoolVar(variable_name);
        }

        for (int k = 1; k < n; k++)
            solver->MakeRowConstraint(LinearExpr(u[k]) >= LinearExpr(u[k + 1])); // (9)

        // (11) obs.
        solver->MakeRowConstraint(LinearExpr(u[0]) == 1.0);
        solver->MakeRowConstraint(LinearExpr(u[n + 1]) == 0.0);

        // flow is "ye" variables for cost arcs
        vector<map<Arc, MPVariable *>> flow(stringpool.size());

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
            {
                string var_name = "flow" + to_string(i) + "_" +
                                  to_string(get<0>(get<0>(arc))) + "_" + to_string(get<1>(get<0>(arc))) + "_" + to_string(get<0>(get<1>(arc))) + "_" + to_string(get<1>(get<1>(arc)));

                internal_flow[arc] = solver->MakeNumVar(0.0, solver->infinity(), var_name);
            }

            // zero_cost_flow is "ye" variables for free slant arcs
            map<FreeSlantArc, MPVariable *> zero_cost_flow;
            for (const auto &arc : free_slant_arcs)
            {
                string var_name = "ye" + to_string(i) + "_" +
                                  to_string(get<0>(get<0>(arc))) + "_" + to_string(get<1>(get<0>(arc))) + "_" + to_string(get<0>(get<1>(arc))) + "_" + to_string(get<1>(get<1>(arc)));

                zero_cost_flow[arc] = solver->MakeNumVar(0.0, solver->infinity(), var_name);
            }

            for (const auto &arc : arcs)
            {
                int j = get<1>(get<1>(arc));

                if (0 < j && j <= n)
                    solver->MakeRowConstraint(LinearExpr(internal_flow[arc]) <= LinearExpr(u[j])); // (10)
            }

            for (const auto &arc : free_slant_arcs)
            {
                int j = get<1>(get<1>(arc));

                if (0 < j && j <= n)
                    solver->MakeRowConstraint(LinearExpr(zero_cost_flow[arc]) <= LinearExpr(u[j])); // (10)
            }

            for (const auto &zcf : zero_cost_flow)
            {
                if (get<2>(zcf.first) == 0)
                {
                    // (5)
                    solver->MakeRowConstraint(
                        LinearExpr(zcf.second) <=
                        1 - LinearExpr(x[get<1>(get<1>(zcf.first)) - 1]));
                }
                else
                {
                    // (6)
                    solver->MakeRowConstraint(
                        LinearExpr(zcf.second) <=
                        LinearExpr(x[get<1>(get<1>(zcf.first)) - 1]));
                }
            }

            LinearExpr source;

            for (const auto &arc : arcs)
                if (get<0>(arc) == make_tuple(0, 0))
                    source += internal_flow[arc];

            for (const auto &arc : free_slant_arcs)
                if (get<0>(arc) == make_tuple(0, 0))
                    source += zero_cost_flow[arc];

            // (3)
            solver->MakeRowConstraint(source == 1.0);

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

                // (4)
                solver->MakeRowConstraint(outflow == inflow);
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

                // (11)
                solver->MakeRowConstraint(inflow - outflow == LinearExpr(u[k]) - LinearExpr(u[k + 1]));
            }

            LinearExpr distance;

            for (const auto &value : internal_flow)
                distance += get<2>(value.first) * LinearExpr(value.second);

            solver->MakeRowConstraint(distance <= d);
        }

        // triangle inequality
        for (int i = 0; i < stringpool.size(); i++)
        {
            LinearExpr distance_i_x;

            for (const auto &value : flow[i])
                distance_i_x += get<2>(value.first) * LinearExpr(value.second);

            for (int j = 0; j < stringpool.size(); j++)
            {
                if (i == j)
                    continue;

                LinearExpr distance_j_x;

                for (const auto &value : flow[j])
                    distance_j_x += get<2>(value.first) * LinearExpr(value.second);

                solver->MakeRowConstraint(dm(i, j) + distance_i_x >= distance_j_x);
            }

            for (int j = i + 1; j < stringpool.size(); j++)
            {
                LinearExpr distance_j_x;

                for (const auto &value : flow[j])
                    distance_j_x += get<2>(value.first) * LinearExpr(value.second);

                solver->MakeRowConstraint(distance_i_x + distance_j_x >= dm(i, j));
            }
        }

        // additional constraints made results worse
        // for (int k = 0; k <= n; k++)
        // {
        //     solver->MakeRowConstraint(higher_str_size - LinearExpr(d) <=
        //                               k - 1 + LinearExpr(u[k]) * (higher_str_size + 1));

        //     solver->MakeRowConstraint(k <=
        //                               p + LinearExpr(d) - higher_str_size * (LinearExpr(u[k]) - 1));
        // }

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
        for (MPVariable *var : u)
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
