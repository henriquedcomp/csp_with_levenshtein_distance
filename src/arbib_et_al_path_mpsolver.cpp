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
    optional<Solution> arbib_et_al_path(
        int n,
        const vector<vector<int>> &stringpool,
        const string_view &a,
        int p)
    {
        unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GUROBI"));
        if (!solver)
        {
            LOG(WARNING) << "Could not create solver GUROBI";
            return nullopt;
        }

        solver->set_time_limit(TIME_LIMIT_IN_MS);

        const double infinity = solver->infinity();
        MPVariable *d = solver->MakeNumVar(0.0, infinity, "d");

        vector<MPVariable *> x;
        solver->MakeNumVarArray(n, 0, a.size() - 1, "x", &x);

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

        for (int i = 0; i < stringpool.size(); i++)
        {
            int nk = stringpool[i].size();

            vector<Arc> arcs;
            GraphConstruction::generate_arcs(nk, n, arcs);

            vector<Vertex> nodes;
            GraphConstruction::generate_nodes(arcs, nodes);

            vector<FreeSlantArc> free_slant_arcs;
            GraphConstruction::generate_free_slant_arcs(stringpool[i], nk, n, free_slant_arcs);

            // flow is "ye" variables for cost arcs
            map<Arc, MPVariable *> flow;
            for (const auto &arc : arcs)
            {
                string var_name = "flow" + to_string(i) + "_" +
                                  to_string(get<0>(get<0>(arc))) + "_" + to_string(get<1>(get<0>(arc))) + "_" + to_string(get<0>(get<1>(arc))) + "_" + to_string(get<1>(get<1>(arc)));

                flow[arc] = solver->MakeNumVar(0.0, solver->infinity(), var_name);
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
                    solver->MakeRowConstraint(LinearExpr(flow[arc]) <= LinearExpr(u[j])); // (10)
            }

            for (const auto &arc : free_slant_arcs)
            {
                int j = get<1>(get<1>(arc));

                if (0 < j && j <= n)
                    solver->MakeRowConstraint(LinearExpr(zero_cost_flow[arc]) <= LinearExpr(u[j])); // (10)
            }

            for (const auto &zcf : zero_cost_flow)
            {
                string name = format("e_{}{}_{}{}", get<0>(get<0>(zcf.first)), get<0>(get<1>(zcf.first)), get<1>(get<0>(zcf.first)), get<1>(get<1>(zcf.first)));

                MPVariable *e = solver->MakeBoolVar(name);

                solver->MakeRowConstraint(get<2>(zcf.first) - LinearExpr(x[get<1>(get<1>(zcf.first)) - 1]) <=
                                          a.size() * LinearExpr(e));

                solver->MakeRowConstraint(LinearExpr(x[get<1>(get<1>(zcf.first)) - 1]) - get<2>(zcf.first) <=
                                          a.size() * LinearExpr(e));

                solver->MakeRowConstraint(LinearExpr(e) + LinearExpr(zcf.second) <= 1);
            }

            LinearExpr source;

            for (const auto &arc : arcs)
                if (get<0>(arc) == make_tuple(0, 0))
                    source += flow[arc];

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

                // (11)
                solver->MakeRowConstraint(inflow - outflow == LinearExpr(u[k]) - LinearExpr(u[k + 1]));
            }

            LinearExpr distance;

            for (const auto &value : flow)
                distance += get<2>(value.first) * LinearExpr(value.second);

            solver->MakeRowConstraint(distance <= d);
        }

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
