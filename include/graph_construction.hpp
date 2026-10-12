#pragma once

#include <vector>
#include <tuple>
#include <set>

#include "config.hpp"

using namespace std;

using Vertex = tuple<int, int>;
using Arc = tuple<Vertex, Vertex, double>;
using FreeSlantArc = tuple<Vertex, Vertex, int>;

class GraphConstruction
{
public:
    static void generate_arcs(int nk, int n, vector<Arc> &arcs)
    {
        for (int i = 0; i <= nk; i++)
        {
            for (int j = 0; j <= n; j++)
            {
                if (i < nk)
                {
                    arcs.push_back(
                        Arc(
                            Vertex(i, j),
                            Vertex(i + 1, j),
                            CDEL));
                }

                if (j < n)
                {
                    arcs.push_back(
                        Arc(
                            Vertex(i, j),
                            Vertex(i, j + 1),
                            CINS));
                }

                if (i < nk && j < n)
                {
                    arcs.push_back(
                        Arc(
                            Vertex(i, j),
                            Vertex(i + 1, j + 1),
                            CSUB));
                }
            }
        }
    }

    static void generate_free_slant_arcs(
        const vector<int> &word,
        int nk, int n,
        vector<FreeSlantArc> &free_slant_arcs)
    {
        for (int i = 0; i < nk; i++)
        {
            for (int j = 0; j < n; j++)
            {
                free_slant_arcs.push_back(
                    FreeSlantArc(
                        Vertex(i, j),
                        Vertex(i + 1, j + 1),
                        word[i]));
            }
        }
    }

    static void generate_nodes(
        vector<Arc> &arcs,
        vector<Vertex> &nodes_list)
    {
        set<Vertex> nodes;

        for (const Arc &arc : arcs)
        {
            if (!nodes.contains(get<0>(arc)))
                nodes.insert(get<0>(arc));

            if (!nodes.contains(get<1>(arc)))
                nodes.insert(get<1>(arc));
        }

        nodes_list.clear();
        nodes_list.reserve(nodes.size());
        nodes_list.assign(nodes.begin(), nodes.end());
    }
};