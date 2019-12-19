#include "movingcell.h"

namespace Common
{
    void MovingCell::reserve_for_mcp(size_t size)
    {
        mcp_.reserve(size);
    }

    void MovingCell::resize_cell()
    {
        if (mcp_.empty()) {
            return;
        }

        /*int limit = 1e6;

        int count = mcp_.size();
        for (int i=0; i<mcp_.size(); ++i)
        {
            int size = mcp_[i].cell_.size();
            int c = std::ceil(size / limit);
            count += c;
        }

        cell_.resize(count);

        for (int i=0; i<count;)
        {
            int size = mcp_[i].cell_.size();
            int c = std::ceil(size / limit);
            int current = size;

            for (int j=0; j<c; ++j)
            {
                int m = std::ceil(current / (c-j));
                cell_[i+j].reserve(m);
                current -= m;
            }

            i += c;
        }*/

        cell_.resize(mcp_.size());
        assert(cell_.size() == mcp_.size());
        int count = 0;
        for (int i=0; i<mcp_.size(); ++i)
        {
            int size = mcp_[i].cell_.size();
            cell_[i].reserve(size);
            count += size;
        }
    }

    const std::map<unsigned int, unsigned int>& MovingCell::index_to_binindex() const
    {
        return index_to_binindex_;
    }

    void MovingCell::add_cell_real(const std::deque<Mesh>& mesh, std::string s)
    {
        if (cell_.empty())
        {
            resize_cell();
        }

        assert(cell_.size() == mcp_.size());

        for (int i=0; i<mcp_.size(); ++i)
        {
            for (int j=0; j<mcp_[i].cell_.size(); ++j)
            {
                /*auto m = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& tm){return tm.tag()() == mcp_[i].mesh_[j];});
                assert(m != mesh.end());
                std::vector<MeshCell>::const_iterator c;
                MeshCell tmc;
                tmc.set_tag(mcp_[i].cell_[j]);
                if (s == "interior")
                {
                    //c = std::find_if(m->cell().begin(), m->cell().end(), [&](const MeshCell& tmc){return tmc.tag()() == mcp_[i].cell_[j];});
                    c = std::lower_bound(m->cell().begin(), m->cell().end(), tmc, [&](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
                    assert(c != m->cell().end());
                    assert(c->tag() == tmc.tag());
                }
                else if (s == "wall")
                {
                    //c = std::find_if(m->wall_boundaries().begin(), m->wall_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag()() == mcp_[i].cell_[j];});
                    c = std::lower_bound(m->wall_boundaries().begin(), m->wall_boundaries().end(), tmc, [&](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
                    assert(c != m->wall_boundaries().end());
                    assert(c->tag() == tmc.tag());
                }
                else if (s == "outer")
                {
                    //c = std::find_if(m->outer_boundaries().begin(), m->outer_boundaries().end(), [&](const MeshCell& tmc){return tmc.tag()() == mcp_[i].cell_[j];});
                    c = std::lower_bound(m->outer_boundaries().begin(), m->outer_boundaries().end(), tmc, [&](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
                    assert(c != m->outer_boundaries().end());
                    assert(c->tag() == tmc.tag());
                }
                else
                {
                    assert(false);
                }*/

                assert(i < cell_.size());
                //cell_[i].push_back(*c);
                assert(mcp_[i].address_[j] != nullptr);
                cell_[i].push_back(*(mcp_[i].address_[j]));

                //auto itt = std::find_if(indices_.begin(), indices_.end(), [&](const std::pair<Tag, Tag>& _dc) { return (_dc.first == c->parent_mesh() && _dc.second == c->tag()); });
                //if (itt == indices_.end()) {
                    //indices_.push_back(std::pair<Tag, Tag>(c->parent_mesh(), c->tag()));
                //}
            }
        }
    }

    void MovingCell::add_cell_pre(unsigned int rank, const MeshCell& c, unsigned int binindex)
    {
        auto it = std::find_if(mcp_.begin(), mcp_.end(), [&](const auto& tmcp){return tmcp.bin_ == binindex;});
        if (it == mcp_.end())
        {
            rank_to_index_.insert(std::pair<unsigned int, unsigned int> (rank, mcp_.size()));
            index_to_rank_.insert(std::pair<unsigned int, unsigned int> (mcp_.size(), rank));
            binindex_to_index_.insert(std::pair<unsigned int, unsigned int> (binindex, mcp_.size()));
            index_to_binindex_.insert(std::pair<unsigned int, unsigned int> (mcp_.size(), binindex));
            mcp_.push_back(MovingCellPre(c.tag()(), c.parent_mesh()(), rank, binindex, &c));
        }
        else {
            it->add(c.tag()(), c.parent_mesh()(), rank, &c);
        }
    }

    void MovingCell::add_cell(unsigned int rank, const MeshCell& c, unsigned int binindex)
    {
        //auto it = rank_to_index_.find(rank);
        auto it = binindex_to_index_.find(binindex);
        //if (it == rank_to_index_.end())
        if (it == binindex_to_index_.end())
        {
            std::vector<MeshCell> temp;
            temp.push_back(c);
            rank_to_index_.insert(std::pair<unsigned int, unsigned int> (rank, cell_.size()));
            index_to_rank_.insert(std::pair<unsigned int, unsigned int> (cell_.size(), rank));
            binindex_to_index_.insert(std::pair<unsigned int, unsigned int> (binindex, cell_.size()));
            index_to_binindex_.insert(std::pair<unsigned int, unsigned int> (cell_.size(), binindex));
            cell_.push_back(temp);
        }
        else
            cell_[it->second].push_back(c);
    }

    unsigned int MovingCell::index_to_rank(unsigned int i)
    {
        auto it = index_to_rank_.find(i);
        assert(it != index_to_rank_.end());
        return it->second;
    }

    unsigned int MovingCell::index_to_binindex(unsigned int i)
    {
        auto it = index_to_binindex_.find(i);
        assert(it != index_to_binindex_.end());
        return it->second;
    }

    const std::vector<std::vector<MeshCell>>& MovingCell::cell() const
    {
        return cell_;
    }

    const std::vector<MeshCell>& MovingCell::cell(unsigned int i) const
    {
        return cell_[i];
    }

    void MovingCell::clear()
    {
        std::vector<std::vector<MeshCell>>().swap(cell_);
        std::vector<MovingCellPre>().swap(mcp_);
        rank_to_index_.clear();
        index_to_rank_.clear();
        binindex_to_index_.clear();
        index_to_binindex_.clear();
    }

    MovingCellPre::MovingCellPre(int cell, int mesh, unsigned int dest, unsigned int bin, const MeshCell* addr): bin_(bin)
    {
        add(cell, mesh, dest, addr);
    }

    void MovingCellPre::add(int cell, int mesh, unsigned int dest, const MeshCell* addr)
    {
        cell_.push_back(cell);
        mesh_.push_back(mesh);
        dest_.push_back(dest);
        address_.push_back(addr);
    }
}
