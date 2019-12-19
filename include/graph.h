#ifndef GRAPH_H
#define GRAPH_H

#include <map>
#include <vector>
#include <memory>
#include <set>
#include <cassert>
#include <algorithm>
#include <iostream>
#include "metis.h"
#include "tag.h"

namespace Common
{
    class Graph
    {
        public:

            Graph(int w0, int w1, int w2, int w3);
            void refine(int vertexid, int quadid, int newquadid, int newvertexid, int w0, int w1, int w2, int w3);
            void connect();
            void print();
            bool partition(idx_t nparts, double& dev, std::vector<std::vector<BinRMTag>>& aug_bin_tag);

        private:

            enum VertexPos
            {
                LL=0,
                LR=1,
                UL=2,
                UR=3
            };

            enum Pos
            {
                L,
                B,
                R,
                T
            };

            struct Quad;

            struct Vertex
            {
                int id_;
                int parent_quad_;
                int pos_;
                int weight_;
                std::unique_ptr<Quad> quad_;
                std::vector<Vertex*> nei_;
            };

            struct Quad
            {
                int id_;
                //std::array<Vertex, 4> vertex_;
                Vertex vertex_[4];

                Quad(int id, int newvertexid, int w0, int w1, int w2, int w3);
            };

            Quad quad_;
            std::map<int, Quad*> id_to_quad_;
            int nedge_;
            int current_max_id_;
            std::vector<Vertex*> vertex_;

            void connect(Vertex& a, Vertex& b, Pos pos);
            void connect(Quad* quad);
            Vertex& vertex(int vertexid, int quadid); // vertexpos such as LL.
            void reset_connection(Quad& quad);
            void reset_connection();
    };
}

#endif
