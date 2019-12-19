#include "mesh_transfer_info.h"

namespace Common
{
    int MeshTransferInfo::Dest::ncell() const
    {
        int count = 0;
        for (const auto& m: mesh_)
        {
            count += m.cell_.size();
        }

        return count;
    }

    int pairing_function(int a, int b)
    {
        int c = (a + b) * (a + b + 1) / 2 + b;

        return c;
    }

    std::vector<MeshCell> MeshTransferInfo::Dest::make_cells(const std::deque<Mesh>& m, const AABB& aabb, bool check_residency) const
    {
        std::vector<MeshCell> mb;

        for (const auto& mm: mesh_)
        {
            const Tag& mt = mm.tag_;
            {
                for (auto itt = mm.cell_.begin(); itt != mm.cell_.end(); ++itt)
                {
                    auto mit = std::find_if(m.begin(), m.end(), [&] (const Mesh& mmm) {return mmm.tag() == mt;});
                    assert(mit != m.end());

                    if (check_residency)
                    {
                        //if (!aabb.do_intersect(mit->cell(*itt).polygon().centroid())) {
                        if (!aabb.do_intersect(mit->cell(*itt).poly().centroid())) {
                            continue;
                        }
                    }

                    mb.push_back(mit->cell(*itt));
                }
            }
        }

        return mb;
    }

    std::deque<Mesh> MeshTransferInfo::Dest::make_meshes_without_pts(const std::deque<Mesh>& m) const
    {
        std::deque<Mesh> mb;

        for (const auto& mm: mesh_)
        {
            const Tag& mt = mm.tag_;
            auto mit = std::find_if(m.begin(), m.end(), [&] (const Mesh& mmm) {return mmm.tag() == mt;});
            assert(mit != m.end());
            auto it = std::find_if(mb.begin(), mb.end(), [&] (const Mesh& mmm) {return mmm.tag() == mt;});
            if (it == mb.end())
            {
                mb.push_back(Mesh());
                Mesh& mbr = mb.back();
                mbr.set_parent_mesh(mt);
                mbr.set_tag(mt);

                // possibly capacity > size.
                mbr.reserve_interior_only(mm.cell_.size());
                mbr.reserve_wall(mit->wall_boundaries().size());
                mbr.reserve_dirichlet(mit->dirichlet_boundaries().size());
                mbr.reserve_empty(mit->empty_boundaries().size());
                mbr.reserve_farfield(mit->farfield_boundaries().size());

                for (auto itt = mm.cell_.begin(); itt != mm.cell_.end(); ++itt)
                {
                    mbr.add_cell_nonsorted(mit->cell(*itt));

                    for (const Tag& t: mit->cell(*itt).wall_boundary())
                    {
                        mbr.add_element_nonsorted(mit->wall_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).dirichlet_boundary())
                    {
                        mbr.add_element_nonsorted(mit->dirichlet_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).empty_boundary())
                    {
                        mbr.add_element_nonsorted(mit->empty_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).farfield_boundary())
                    {
                        mbr.add_element_nonsorted(mit->farfield_boundary(t));
                    }
                }

                mbr.shrink_cells();
                mbr.sort_cells();
                //mbr.set_cell_tag_index_map();
            }
        }
        
        return mb;
    }

    std::deque<Mesh> MeshTransferInfo::Dest::make_meshes(const std::deque<Mesh>& m, const AABB& aabb, bool check_residency) const
    {
        std::deque<Mesh> mb;

        for (const auto& mm: mesh_)
        {
            const Tag& mt = mm.tag_;
            auto mit = std::find_if(m.begin(), m.end(), [&] (const Mesh& mmm) {return mmm.tag() == mt;});
            assert(mit != m.end());
            if (mit->wall_boundaries().size() == 0)
            {
                for (const MeshCell& mc: mit->cell())
                {
                    assert(mc.wall_boundary().empty());
                }
            }
            auto it = std::find_if(mb.begin(), mb.end(), [&] (const Mesh& mmm) {return mmm.tag() == mt;});
            if (it == mb.end())
            {
                mb.push_back(Mesh());
                Mesh& mbr = mb.back();
                mbr.set_parent_mesh(mt);
                mbr.set_tag(mt);

                // possibly capacity > size.
                mbr.reserve_interior(mm.cell_.size());
                mbr.reserve_wall(mit->wall_boundaries().size()); // notice argument difference between upper line.
                mbr.reserve_dirichlet(mit->dirichlet_boundaries().size());
                mbr.reserve_empty(mit->empty_boundaries().size());
                mbr.reserve_farfield(mit->farfield_boundaries().size());

                for (auto itt = mm.cell_.begin(); itt != mm.cell_.end(); ++itt)
                {
                    if (check_residency)
                    {
                        if (!aabb.do_intersect(mit->cell(*itt).poly().centroid())) {
                            continue;
                        }
                    }

                    mbr.add_interior_nonsorted_addpoint(mit->cell(*itt));

                    for (const Tag& t: mit->cell(*itt).wall_boundary())
                    {
                        assert(mit->wall_boundaries().size() != 0);
                        mbr.add_element(mit->wall_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).dirichlet_boundary())
                    {
                        assert(mit->dirichlet_boundaries().size() != 0);
                        mbr.add_element(mit->dirichlet_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).empty_boundary())
                    {
                        assert(mit->empty_boundaries().size() != 0);
                        mbr.add_element(mit->empty_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).farfield_boundary())
                    {
                        assert(mit->farfield_boundaries().size() != 0);
                        mbr.add_element(mit->farfield_boundary(t));
                    }
                }

                mbr.remove_duplicate_cells();
                mbr.remove_duplicate_points();
                mbr.shrink_cells();
                mbr.shrink_points();
                mbr.sort_cells();
                mbr.sort_points();
                mbr.remove_parent_cells_of_vertices_of_all_cells();
                mbr.set_parent_cell_of_vertices_of_all_cells();
                mbr.remove_parent_cells_of_all_points();
                mbr.update_points_from_cell_vertices();

                //mbr.set_cell_tag_index_map();
                //mbr.set_point_tag_index_map();
            }
            else
            {
                assert(false);
                for (auto itt = mm.cell_.begin(); itt != mm.cell_.end(); ++itt)
                {
                    if (check_residency)
                    {
                        if (!aabb.do_intersect(mit->cell(*itt).poly().centroid())) {
                            continue;
                        }
                    }

                    assert(!it->query(*itt));
                    assert(mit->cell(*itt).tag().isvalid());

                    //it->add_element_no_check(mit->cell(*itt));
                    it->add_element_nonsorted(mit->cell(*itt));

                    for (const Tag& t: mit->cell(*itt).wall_boundary())
                    {
                        assert(mit->wall_boundaries().size() != 0);
                        it->add_element_nonsorted(mit->wall_boundary(t));
                        //it->add_element_no_check(mit->wall_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).dirichlet_boundary())
                    {
                        assert(mit->dirichlet_boundaries().size() != 0);
                        it->add_element_nonsorted(mit->dirichlet_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).empty_boundary())
                    {
                        assert(mit->empty_boundaries().size() != 0);
                        it->add_element_nonsorted(mit->empty_boundary(t));
                    }
                    for (const Tag& t: mit->cell(*itt).farfield_boundary())
                    {
                        assert(mit->farfield_boundaries().size() != 0);
                        it->add_element_nonsorted(mit->farfield_boundary(t));
                    }
                }
            }
        }
        
        for (Mesh& m: mb)
        {
            //m.remove_dup_cells_and_points();
            m.shrink_points();
            m.sort_cells();
            //m.set_cell_tag_index_map();
            //m.set_point_tag_index_map();
        }

        return mb;
    }

    const std::vector<int> MeshTransferInfo::rank() const
    {
        return rank_;
    }

    const std::vector<MeshTransferInfo::Dest>& MeshTransferInfo::dest() const
    {
        return dest_;
    }

    MeshTransferInfo::MTIMesh::MTIMesh(const Tag& tag, const Tag& ctag): tag_(tag)
    {
        //cell_.push_front(ctag);
        cell_.push_back(ctag);
    }

    MeshTransferInfo::Dest::Dest(int rank, const Tag& mtag, const Tag& ctag, const BinRMTag& sp_tag): rank_(rank), sp_tag_(sp_tag)
    {
        assert(sp_tag.isvalid());
        auto it = std::find_if(mesh_.begin(), mesh_.end(), [&] (const MTIMesh& m) {return m.tag_ == mtag;});
        if (it == mesh_.end())
        {
            //mesh_.push_front(MTIMesh(mtag, ctag));
            mesh_.push_back(MTIMesh(mtag, ctag));
        }
        else
        {
            //it->cell_.push_front(ctag);
            it->cell_.push_back(ctag);
        }
    }

    MeshTransferInfo::MeshTransferInfo(int source, int dest, const Tag& mtag, const Tag& ctag, const BinRMTag& sp_tag): source_(source)
    {
        assert(sp_tag.isvalid());
        add(dest, mtag, ctag, sp_tag);
    }

    MeshTransferInfo::MeshTransferInfo(int source): source_(source)
    {
    }

    void MeshTransferInfo::Dest::add(const Tag& mtag, const Tag& ctag)
    {
        auto it = std::find_if(mesh_.begin(), mesh_.end(), [&] (const MTIMesh& m) {return m.tag_ == mtag;});
        if (it == mesh_.end())
        {
            mesh_.push_back(MTIMesh(mtag, ctag));
        }
        else
        {
            it->cell_.push_back(ctag);
        }
    }

    bool MeshTransferInfo::bin_exist(const BinRMTag& sp_tag) const
    {
        assert(sp_tag.isvalid());
        int c = std::count_if(dest_.begin(), dest_.end(), [&] (const Dest& dst) {return dst.sp_tag_ == sp_tag;});

        if (c == 0)
        {
            return false;
        }

        return true;
    }

    bool MeshTransferInfo::dest_exist(int rank) const
    {
        int c = std::count_if(dest_.begin(), dest_.end(), [&] (const Dest& dst) {return dst.rank_ == rank;});

        if (c == 0)
        {
            return false;
        }

        return true;
    }

    void MeshTransferInfo::add(int dest, const Tag& mtag, const Tag& ctag, const BinRMTag& sp_tag)
    {
        assert(sp_tag.isvalid());
        int c = std::count(rank_.begin(), rank_.end(), dest);
        if (c == 0)
        {
            rank_.push_back(dest);
        }

        auto it = std::find_if(dest_.begin(), dest_.end(), [&] (const Dest& dst) {return dst.sp_tag_ == sp_tag;});
        if (it == dest_.end())
        {
            dest_.push_back(Dest(dest, mtag, ctag, sp_tag));
        }
        else
        {
            if (it->rank_ != dest)
            {
                std::cout << "rank: " << it->rank_ << std::endl;
                std::cout << "dest: " << dest << std::endl;
                std::cout << "sp_tag.rmtag(): " << sp_tag.rmtag()() << std::endl;
                std::cout << "sp_tag.bintag(): " << sp_tag.bintag()() << std::endl;
            }
            assert(it->rank_ == dest);
            it->add(mtag, ctag);
        }
    }

    int MeshTransferInfo::source() const
    {
        return source_;
    }
}
