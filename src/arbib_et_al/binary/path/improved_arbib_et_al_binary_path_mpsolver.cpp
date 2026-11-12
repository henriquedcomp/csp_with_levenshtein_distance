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
        int sol_size_ub,
        int sol_size_lb,
        const vector<vector<int>> &input_strings,
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
        solver->MakeBoolVarArray(sol_size_ub, "t", &t);

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

        // flow is "ye" variables for cost arcs
        vector<map<Arc, MPVariable *>> flow(input_strings.size());

        for (int i = 0; i < input_strings.size(); i++)
        {
            int nk = input_strings[i].size();

            vector<Arc> arcs;
            GraphConstruction::generate_arcs(nk, sol_size_ub, arcs);

            vector<Vertex> nodes;
            GraphConstruction::generate_nodes(arcs, nodes);

            vector<FreeSlantArc> free_slant_arcs;
            GraphConstruction::generate_free_slant_arcs(input_strings[i], nk, sol_size_ub, free_slant_arcs);

            auto &internal_flow = flow[i];
            for (const auto &arc : arcs)
            {
                string var_name = "flow" + to_string(i) + "_" +
                                  to_string(get<0>(get<0>(arc))) + "_" + to_string(get<1>(get<0>(arc))) + "_" + to_string(get<0>(get<1>(arc))) + "_" + to_string(get<1>(get<1>(arc)));

                internal_flow[arc] = solver->MakeNumVar(0.0, 1.0, var_name);
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
                    solver->MakeRowConstraint(LinearExpr(internal_flow[arc]) <= LinearExpr(l[j])); // (4.39)
            }

            for (const auto &arc : free_slant_arcs)
            {
                int j = get<1>(get<1>(arc));

                if (sol_size_lb < j && j <= sol_size_ub)
                    solver->MakeRowConstraint(LinearExpr(zero_cost_flow[arc]) <= LinearExpr(l[j])); // (4.39)
            }

            for (const auto &zcf : zero_cost_flow)
            {
                if (get<2>(zcf.first) == 0)
                {
                    solver->MakeRowConstraint(
                        LinearExpr(zcf.second) <=
                        1 - LinearExpr(t[get<1>(get<1>(zcf.first)) - 1])); // (4.34)
                }
                else
                {
                    solver->MakeRowConstraint(
                        LinearExpr(zcf.second) <=
                        LinearExpr(t[get<1>(get<1>(zcf.first)) - 1])); // (4.33)
                }
            }

            LinearExpr source;

            for (const auto &arc : arcs)
                if (get<0>(arc) == make_tuple(0, 0))
                    source += internal_flow[arc];

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

                solver->MakeRowConstraint(outflow == inflow); // (4.32)
            }

            for (int k = sol_size_lb; k <= sol_size_ub; k++)
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

                solver->MakeRowConstraint(inflow - outflow == LinearExpr(l[k]) - LinearExpr(l[k + 1])); // (4.38)
            }

            LinearExpr distance;

            for (const auto &value : internal_flow)
                distance += get<2>(value.first) * LinearExpr(value.second);

            solver->MakeRowConstraint(distance <= r); // (4.30)
        }

        // additional constraints
        LinearExpr sum_u = sol_size_lb;
        for (int k = sol_size_lb + 1; k <= sol_size_ub; k++)
            sum_u += LinearExpr(l[k]);

        solver->MakeRowConstraint(CINS * higher_str_size - LinearExpr(r) <=
                                  CINS * sum_u); // (4.43)

        solver->MakeRowConstraint(CDEL * sum_u <=
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
