#include "graph.h"

namespace Common
{
    Graph::Vertex& Graph::vertex(int vertexid, int quadid)
    {
        auto it = id_to_quad_.find(quadid);
        assert(it != id_to_quad_.end());
        Quad* quadit = it->second;
        assert(quadit != nullptr);
        for (int i=0; i<4; ++i)
        {
            if (quadit->vertex_[i].id_ == vertexid)
            {
                return quadit->vertex_[i];
            }
        }
        //auto vertexit = std::find_if(quadit->vertex_.begin(), quadit->vertex_.begin(), [&](const Vertex& v){return vertexpos == v.;});
        //assert(vertexit != quadit->vertex_.end());
        //return quadit->vertex_[vertexpos];
    }

    //Graph::Quad::Quad(int id, int newvertexid, int w0, int w1, int w2, int w3): id_(id)
    Graph::Quad::Quad(int id, int newvertexid, int w0, int w1, int w2, int w3): id_(id)
    {
        vertex_[LL].parent_quad_ = id_;
        vertex_[LR].parent_quad_ = id_;
        vertex_[UR].parent_quad_ = id_;
        vertex_[UL].parent_quad_ = id_;

        //vertex_[LL].id_ = LL;
        //vertex_[LR].id_ = LR;
        //vertex_[UR].id_ = UR;
        //vertex_[UL].id_ = UL;

        vertex_[LL].id_ = newvertexid;
        vertex_[LR].id_ = newvertexid + 1;
        vertex_[UL].id_ = newvertexid + 2;
        vertex_[UR].id_ = newvertexid + 3;

        vertex_[LL].pos_ = LL;
        vertex_[LR].pos_ = LR;
        vertex_[UR].pos_ = UR;
        vertex_[UL].pos_ = UL;

        vertex_[LL].weight_ = w0;
        vertex_[LR].weight_ = w1;
        vertex_[UL].weight_ = w2;
        vertex_[UR].weight_ = w3;

        //vertex_[LL].area_ = a0;
        //vertex_[LR].area_ = a1;
        //vertex_[UL].area_ = a2;
        //vertex_[UR].area_ = a3;
    }

    void Graph::reset_connection(Quad& quad)
    {
        for (int i=0; i<4; ++i)
        {
            quad.vertex_[i].nei_.clear();
            //for (Vertex*& v: quad.vertex_[i].nei_) {
                //v = nullptr;
            //}
        }

        for (int i=0; i<4; ++i)
        {
            if (quad.vertex_[i].quad_) {
                reset_connection(*(quad.vertex_[i].quad_.get()));
            }
        }
    }

    void Graph::reset_connection()
    {
        reset_connection(quad_);
    }

    Graph::Graph(int w0, int w1, int w2, int w3): current_max_id_(0), quad_(0,0, w0, w1, w2, w3), nedge_(0)
    {
        // weights are assigned in this order: LL, LR, UL, UR

        id_to_quad_.insert(std::make_pair(0, &quad_));
    }

    void Graph::connect(Quad* quad)
    {
        if (quad == nullptr) {
            return;
        }

        if (quad->vertex_[LL].quad_) {
            connect(quad->vertex_[LL].quad_.get());
        }
        if (quad->vertex_[LR].quad_) {
            connect(quad->vertex_[LR].quad_.get());
        }
        if (quad->vertex_[UL].quad_) {
            connect(quad->vertex_[UL].quad_.get());
        }
        if (quad->vertex_[UR].quad_) {
            connect(quad->vertex_[UR].quad_.get());
        }

        connect(quad->vertex_[LL], quad->vertex_[LR], L);
        connect(quad->vertex_[LR], quad->vertex_[UR], B);
        connect(quad->vertex_[UR], quad->vertex_[UL], R);
        connect(quad->vertex_[UL], quad->vertex_[LL], T);
    }

    void Graph::connect()
    {
        vertex_.clear();
        connect(&quad_);
    }

    void Graph::connect(Vertex& a, Vertex& b, Pos pos)
    {
        if (!a.quad_)
        {
            if (!b.quad_) {
                a.nei_.push_back(&b);
                b.nei_.push_back(&a);
                ++nedge_;
                //if (a.weight_ != 0.)
                {
                    auto it = std::find_if(vertex_.begin(), vertex_.end(), [&](const Vertex* v){return v->id_ == a.id_;});
                    if (it == vertex_.end()) {
                        vertex_.push_back(&a);
                    }
                }
                //if (b.weight_ != 0.)
                {
                    auto it = std::find_if(vertex_.begin(), vertex_.end(), [&](const Vertex* v){return v->id_ == b.id_;});
                    if (it == vertex_.end()) {
                        vertex_.push_back(&b);
                    }
                }
            }
            else
            {
                if (pos == L)
                {
                    connect(a, b.quad_->vertex_[LL], pos);
                    connect(a, b.quad_->vertex_[UL], pos);
                }
                else if (pos == B)
                {
                    connect(a, b.quad_->vertex_[LL], pos);
                    connect(a, b.quad_->vertex_[LR], pos);
                }
                else if (pos == R)
                {
                    connect(a, b.quad_->vertex_[LR], pos);
                    connect(a, b.quad_->vertex_[UR], pos);
                }
                else if (pos == T)
                {
                    connect(a, b.quad_->vertex_[UL], pos);
                    connect(a, b.quad_->vertex_[UR], pos);
                }
            }
        }
        else
        {
            if (!b.quad_)
            {
                if (pos == L)
                {
                    connect(a.quad_->vertex_[LR], b, pos);
                    connect(a.quad_->vertex_[UR], b, pos);
                }
                else if (pos == B)
                {
                    connect(a.quad_->vertex_[UL], b, pos);
                    connect(a.quad_->vertex_[UR], b, pos);
                }
                else if (pos == R)
                {
                    connect(a.quad_->vertex_[UL], b, pos);
                    connect(a.quad_->vertex_[LL], b, pos);
                }
                else if (pos == T)
                {
                    connect(a.quad_->vertex_[LL], b, pos);
                    connect(a.quad_->vertex_[LR], b, pos);
                }
            }
            else
            {
                if (pos == L)
                {
                    connect(a.quad_->vertex_[LR], b.quad_->vertex_[LL], pos);
                    connect(a.quad_->vertex_[UR], b.quad_->vertex_[UL], pos);
                }
                else if (pos == B)
                {
                    connect(a.quad_->vertex_[UL], b.quad_->vertex_[LL], pos);
                    connect(a.quad_->vertex_[UR], b.quad_->vertex_[LR], pos);
                }
                else if (pos == R)
                {
                    connect(a.quad_->vertex_[LL], b.quad_->vertex_[LR], pos);
                    connect(a.quad_->vertex_[UL], b.quad_->vertex_[UR], pos);
                }
                else if (pos == T)
                {
                    connect(a.quad_->vertex_[LL], b.quad_->vertex_[UL], pos);
                    connect(a.quad_->vertex_[LR], b.quad_->vertex_[UR], pos);
                }
            }
        }
    }

    void Graph::refine(int vertexid, int quadid, int newquadid, int newvertexid, int w0, int w1, int w2, int w3)
    {
        reset_connection();
        Vertex& vtx = vertex(vertexid, quadid);
        assert(vtx.id_ == vertexid);
        
        auto it = id_to_quad_.find(quadid);
        assert(it != id_to_quad_.end());
        auto quadit = it->second;

        assert(vtx.quad_ == nullptr);
        vtx.quad_ = std::make_unique<Quad>(Quad(newquadid, newvertexid, w0, w1, w2, w3));
        //std::cout << "vtx.quad_->id_: " << vtx.quad_->id_ << std::endl;
        //std::cout << "w0: " << w0 << std::endl;
        //std::cout << "w1: " << w1 << std::endl;
        //std::cout << "w2: " << w2 << std::endl;
        //std::cout << "w3: " << w3 << std::endl;
        //std::cout << "weight0: " << vtx.quad_->vertex_[LL].weight_ << std::endl;
        //std::cout << "weight1: " << vtx.quad_->vertex_[LR].weight_ << std::endl;
        //std::cout << "weight2: " << vtx.quad_->vertex_[UL].weight_ << std::endl;
        //std::cout << "weight2: " << vtx.quad_->vertex_[UR].weight_ << std::endl;
        //id_to_quad_.insert(std::make_pair(newquadid->first+1, vtx.quad_.get()));
        id_to_quad_.insert(std::make_pair(newquadid, vtx.quad_.get()));
    }

    void Graph::print()
    {
        std::cout << vertex_.size() << " " << nedge_ << std::endl;
        for (Vertex* v: vertex_)
        {
            std::cout << v->id_;
            for (Vertex* n: v->nei_)
            {
                std::cout << " ";
                std::cout << n->id_;
            }
            std::cout << std::endl;
        }
    }
    
    bool Graph::partition(idx_t nparts, double& dev, std::vector<std::vector<BinRMTag>>& aug_bin_tag)
    {
        bool balance = true;

        idx_t nvtxs = vertex_.size();
        idx_t ncon = 1;
        //idx_t ncon = 2; // multi-constraint: 1) load of bin, 2) area of bin

        idx_t xadj[nvtxs+1];
        idx_t adjncy[nedge_ * 2];
        idx_t vwgt[nvtxs * ncon];
        xadj[0] = 0;
        int i=0;
        int j=0;
        for (Vertex* v: vertex_)
        {
            vwgt[i] = v->weight_;
            //vwgt[i+1] = v->area_;
            ++i;
            //i += ncon;
        }
        i=1;
        for (Vertex* v: vertex_)
        {
            int neisize = 0;
            for (Vertex* n: v->nei_)
            {
                assert(n != nullptr);
                auto it = std::find_if(vertex_.begin(), vertex_.end(), [&](const Vertex* v){return v->id_ == n->id_;});
                if (it == vertex_.end()) {
                    continue;
                }
                adjncy[j] = std::distance(vertex_.begin(), it);
                ++j;
                ++neisize;
            }
            xadj[i] = xadj[i-1] + neisize;
            ++i;
        }

        /*std::cout << "nvtxs: " << nvtxs << std::endl;
        std::cout << "ncon: " << ncon << std::endl;
        std::cout << "nparts: " << nparts << std::endl;

        std::cout << "xadj: ";
        for (int i=0; i<nvtxs+1; ++i) {
            std::cout << xadj[i] << " ";
        }

        std::cout << std::endl;

        std::cout << "adjncy: ";
        for (int i=0; i<nedge_*2; ++i) {
            std::cout << adjncy[i] << " ";
        }

        std::cout << std::endl;

        std::cout << "vwgt: ";
        for (int i=0; i<nvtxs; ++i) {
            std::cout << vwgt[i] << " ";
        }

        std::cout << std::endl;*/

        idx_t obj_val;
        idx_t part[nvtxs];


        //std::cout << "nvertices: " << vertex_.size() << std::endl;
        //std::cout << "nparts: " << nparts << std::endl;

        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_CONTIG] = 1;
        //options[METIS OPTION OBJTYPE] = METIS_OBJTYPE_VOL;

        int res = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, NULL, options, &obj_val, part);
        if (res != METIS_OK) {
            balance = false;
        }

        //if (balance)
        //{
            //for (int i=0; i<nvtxs; ++i) {
                //std::cout << vertex_[i]->parent_quad_ << " " << vertex_[i]->id_ << " " << vertex_[i]->weight_ << " " << part[i] << std::endl;
            //}
        //}

        std::vector<int> group(nparts, 0);

        for (int i=0; i<nvtxs; ++i) {
            group[part[i]] += vertex_[i]->weight_;
        }
        //for (int i=0; i<nparts; ++i) {
            //std::cout << "group[i]: " << group[i] << std::endl;
        //}

        auto result = std::minmax_element(group.begin(), group.end());
        double minload = *result.first;
        double maxload = *result.second;
        if (maxload == 0)
        {
            for (int i=0; i<nvtxs; ++i) {
                std::cout << vertex_[i]->parent_quad_ << " " << vertex_[i]->id_ << " " << vertex_[i]->weight_ << " " << part[i] << std::endl;
            }
        }
        assert(maxload != 0);
        dev = ((maxload - minload) / maxload) * 100.;
        //std::cout << "min: " << minload << std::endl;
        //std::cout << "max: " << maxload << std::endl;
        //std::cout << "dev: " << dev << std::endl;
        //std::cout << "balance: " << balance << std::endl;

        aug_bin_tag.clear();
        aug_bin_tag.resize(nparts);
        for (int i=0; i<nvtxs; ++i)
        {
            //std::cout << "rm: " <<  vertex_[i]->parent_quad_ << std::endl;
            //std::cout << "bt: " <<  vertex_[i]->id_<< std::endl;
            BinRMTag tag(Tag(vertex_[i]->id_), Tag(vertex_[i]->parent_quad_));
            aug_bin_tag[part[i]].push_back(tag);
        }
        assert(!aug_bin_tag.empty());

        /*if (balance)
        {
            std::cout << "aug size: " << aug_bin_tag.size() << std::endl;
            for (int t=0; t<aug_bin_tag.size(); ++t)
            {
                std::cout << "t: " << t << std::endl;
                for (const auto& tt: aug_bin_tag[t])
                {
                    std::cout << "aug[" << t << "].rm: " << tt.rmtag()() << std::endl;
                    std::cout << "aug[" << t << "].bt: " << tt.bintag()() << std::endl;
                }
            }
        }*/

        //std::cout << "done with partitioning" << std::endl;
        return balance;
    }
}
