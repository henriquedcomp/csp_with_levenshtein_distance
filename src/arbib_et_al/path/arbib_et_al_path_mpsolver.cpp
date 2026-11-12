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
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
        const string_view &alphabet,
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

        OffsetVector<MPVariable *> l(sol_size_lb, sol_size_ub + 1);
        for (int k = sol_size_lb; k <= sol_size_ub + 1; k++)
        {
            string variable_name = format("l_{}", k);
            l[k] = solver->MakeBoolVar(variable_name);
        }

        for (int k = sol_size_lb + 1; k < sol_size_ub; k++)
            solver->MakeRowConstraint(LinearExpr(l[k]) >= LinearExpr(l[k + 1])); // (4.35)

        solver->MakeRowConstraint(LinearExpr(l[sol_size_lb]) == 1.0);     // (4.36)
        solver->MakeRowConstraint(LinearExpr(l[sol_size_ub + 1]) == 0.0); // (4.37)

        for (int i = 0; i < input_strings.size(); i++)
        {
            int nk = input_strings[i].size();

            vector<Arc> arcs;
            GraphConstruction::generate_arcs(nk, sol_size_ub, arcs);

            vector<Vertex> nodes;
            GraphConstruction::generate_nodes(arcs, nodes);

            vector<FreeSlantArc> free_slant_arcs;
            GraphConstruction::generate_free_slant_arcs(input_strings[i], nk, sol_size_ub, free_slant_arcs);

            // flow is "ye" variables for cost arcs
            map<Arc, MPVariable *> flow;
            for (const auto &arc : arcs)
            {
                string var_name = "flow" + to_string(i) + "_" +
                                  to_string(get<0>(get<0>(arc))) + "_" + to_string(get<1>(get<0>(arc))) + "_" + to_string(get<0>(get<1>(arc))) + "_" + to_string(get<1>(get<1>(arc)));

                flow[arc] = solver->MakeNumVar(0.0, 1.0, var_name);
            }

            // zero_cost_flow is "ye" variables for free slant arcs
            map<FreeSlantArc, MPVariable *> zero_cost_flow;
            for (const auto &arc : free_slant_arcs)
            {
                string var_name = "ye" + to_string(i) + "_" +
                                  to_string(get<0>(get<0>(arc))) + "_" + to_string(get<1>(get<0>(arc))) + "_" + to_string(get<0>(get<1>(arc))) + "_" + to_string(get<1>(get<1>(arc)));

                zero_cost_flow[arc] = solver->MakeNumVar(0.0, 1.0, var_name);
            }

            for (const auto &arc : arcs)
            {
                int j = get<1>(get<1>(arc));

                if (sol_size_lb < j && j <= sol_size_ub)
                    solver->MakeRowConstraint(LinearExpr(flow[arc]) <= LinearExpr(l[j])); // (4.39)
            }

            for (const auto &arc : free_slant_arcs)
            {
                int j = get<1>(get<1>(arc));

                if (sol_size_lb < j && j <= sol_size_ub)
                    solver->MakeRowConstraint(LinearExpr(zero_cost_flow[arc]) <= LinearExpr(l[j])); // (4.39)
            }

            for (const auto &zcf : zero_cost_flow)
            {
                string name = format("e_{}{}_{}{}", get<0>(get<0>(zcf.first)), get<0>(get<1>(zcf.first)), get<1>(get<0>(zcf.first)), get<1>(get<1>(zcf.first)));

                MPVariable *e = solver->MakeBoolVar(name);

                solver->MakeRowConstraint(get<2>(zcf.first) - LinearExpr(t[get<1>(get<1>(zcf.first)) - 1]) <=
                                          alphabet.size() * LinearExpr(e)); // (4.43)

                solver->MakeRowConstraint(LinearExpr(t[get<1>(get<1>(zcf.first)) - 1]) - get<2>(zcf.first) <=
                                          alphabet.size() * LinearExpr(e)); // (4.44)

                solver->MakeRowConstraint(LinearExpr(e) + LinearExpr(zcf.second) <= 1); // (4.45)
            }

            LinearExpr source;

            for (const auto &arc : arcs)
                if (get<0>(arc) == make_tuple(0, 0))
                    source += flow[arc];

            for (const auto &arc : free_slant_arcs)
                if (get<0>(arc) == make_tuple(0, 0))
                    source += zero_cost_flow[arc];

            solver->MakeRowConstraint(source == 1.0); // (4.31)

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

                solver->MakeRowConstraint(outflow == inflow); // (4.32)
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

                solver->MakeRowConstraint(inflow - outflow == LinearExpr(l[k]) - LinearExpr(l[k + 1])); // (4.38)
            }

            LinearExpr distance;

            for (const auto &value : flow)
                distance += get<2>(value.first) * LinearExpr(value.second);

            solver->MakeRowConstraint(distance <= r); // (4.30)
        }

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
                solution_length = sol_size_lb + i - 1;
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
            absl::ToDoubleSeconds(solver->DurationSinceConstruction()),
            objective->BestBound(),
            (objective->Value() - objective->BestBound()) / objective->Value() * 100.0};

        return solution;
    }
}
