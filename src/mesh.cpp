#include "commonmesh.h"

namespace Common
{
    void Mesh::print_as_vtk(std::string file_name) const
    {
    }

    void Mesh::set_all_cells_as_interior()
    {
        for (MeshCell& mc: cell_)
        {
            mc.set_boutype(boundary_t("interior"));
        }
    }

    void Mesh::print_donor_info() const
    {
        std::ofstream out;
        out.open("donor-info.dat");
        
        for (const MeshCell& mc: cell_)
        {
            out << tag_();
            out << " ";
            out << mc.tag()();
            out << " ";
            out << mc.donor().first();
            out << " ";
            out << mc.donor().second();
        }

        out.close();
    }

    const mfc& Mesh::face() const
    {
        return face_;
    }

    void Mesh::bbox(vec3<double>& min_, vec3<double>& max_) const
    {
        min_.set_x(BIG_POS_NUM);
        min_.set_y(BIG_POS_NUM);
        min_.set_z(BIG_POS_NUM);
        max_.set_x(BIG_NEG_NUM);
        max_.set_y(BIG_NEG_NUM);
        max_.set_z(BIG_NEG_NUM);

        for (const MeshPoint& mp: point_)
        {
            min_.set_x(std::min(min_(0), mp.p().r(0)));
            min_.set_y(std::min(min_(1), mp.p().r(1)));
            min_.set_z(std::min(min_(2), mp.p().r(2)));

            max_.set_x(std::max(max_(0), mp.p().r(0)));
            max_.set_y(std::max(max_(1), mp.p().r(1)));
            max_.set_z(std::max(max_(2), mp.p().r(2)));
        }
    }

    bool Mesh::operator<(const Mesh& other) const
    {
        return tag_ < other.tag();
    }

    size_t Mesh::npartition() const
    {
        // Mesh cells should with set_partition() in connect_cells() before this function.
        
        size_t count = 0;
        for (const MeshCell& mc: cell_)
        {
            if (mc.boutype() == "partition")
            {
                ++count;
            }
        }

        return count;
    }


    /*void Mesh::batch_merge(Mesh& m)
    {
        {
            int size = cell_.size();
            for (auto mc = m.cell().rbegin(); mc != m.cell().rend(); ++mc)
            {
                auto tit = m.cell_tag_index_map.left.find(mc->tag()());
                assert(tit != m.cell_tag_index_map.left.end());
                int rep = tit->second + size;
                bool successful_replace = m.cell_tag_index_map.left.replace_data(tit, rep);
                assert(successful_replace);
            }

            cell_tag_index_map.insert(m.cell_tag_index_map.begin(), m.cell_tag_index_map.end());
        }
        {
            int size = wall_boundaries_.size();
            for (auto mc = m.wall_boundaries().rbegin(); mc != m.wall_boundaries().rend(); ++mc)
            {
                auto tit = m.wall_tag_index_map_.left.find(mc->tag()());
                assert(tit != m.wall_tag_index_map_.left.end());
                int rep = tit->second + size;
                bool successful_replace = m.wall_tag_index_map_.left.replace_data(tit, rep);
                assert(successful_replace);
            }

            wall_tag_index_map_.insert(m.wall_tag_index_map_.begin(), m.wall_tag_index_map_.end());
        }
        {
            int size = outer_boundaries_.size();
            for (auto mc = m.outer_boundaries().rbegin(); mc != m.outer_boundaries().rend(); ++mc)
            {
                auto tit = m.outer_tag_index_map_.left.find(mc->tag()());
                assert(tit != m.outer_tag_index_map_.left.end());
                int rep = tit->second + size;
                bool successful_replace = m.outer_tag_index_map_.left.replace_data(tit, rep);
                assert(successful_replace);
            }

            outer_tag_index_map_.insert(m.outer_tag_index_map_.begin(), m.outer_tag_index_map_.end());
        }

        cell_.insert(cell_.end(), m.cell().begin(), m.cell().end());
        wall_boundaries_.insert(wall_boundaries_.end(), m.wall_boundaries().begin(), m.wall_boundaries().end());
        outer_boundaries_.insert(outer_boundaries_.end(), m.outer_boundaries().begin(), m.outer_boundaries().end());

        int psize = point_.size();
        point_.resize(psize + m.point().size());
        auto it = point_.begin() + psize;
        std::copy_if(m.point().begin(), m.point().end(), it, [&](const MeshPoint& mp){return query_point(mp.tag()) == nullptr;});

        for (auto it = point_.begin() + psize; it != point_.end(); ++it)
        {
            point_tag_index_map.insert(boost::bimap<int, int>::value_type(it->tag()(), std::distance(point_.begin(), it)));
        }

        // now there should not be any duplicates.
        //std::sort(cell_.begin(), cell_.end(), [](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
        //auto itt = std::adjacent_find(cell_.begin(), cell_.end(), [](const MeshCell& left, const MeshCell& right){return left.tag() == right.tag();});
        //assert(itt == cell_.end());
    }*/

    const bimap_int& Mesh::dirichlet_tag_index_map() const
    {
        return dirichlet_tag_index_map_;
    }
        
    void Mesh::destroy_cell_hood()
    {
        for (MeshCell& mc: cell_)
        {
            mc.remove_all_neighbors();
            assert(mc.pnei().empty());
        }
    }

    const std::vector<MeshCell>& Mesh::wall_boundaries() const
    {
        return wall_boundaries_;
    }

    const std::vector<MeshCell>& Mesh::dirichlet_boundaries() const
    {
        return dirichlet_boundaries_;
    }

    const std::vector<MeshCell>& Mesh::empty_boundaries() const
    {
        return wall_boundaries_;
    }

    const std::vector<MeshCell>& Mesh::farfield_boundaries() const
    {
        return wall_boundaries_;
    }

    const MeshCell& Mesh::wall_boundary(const Tag& t) const
    {
        auto it = std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        if(it == wall_boundaries_.end())
        {
            std::cout << "t: " << t() << std::endl;
            std::cout << "size: " << wall_boundaries_.size() << std::endl;
            std::cout << "mesh tag: " << tag_() << std::endl;
        }
        assert(it != wall_boundaries_.end());
        return *it;
    }

    const MeshCell& Mesh::dirichlet_boundary(const Tag& t) const
    {
        auto it = std::find_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        assert(it != dirichlet_boundaries_.end());
        return *it;
    }

    MeshCell& Mesh::wall_boundary_p(const Tag& t)
    {
        auto it = std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        assert(it != wall_boundaries_.end());
        return *it;
    }

    MeshCell& Mesh::dirichlet_boundary_p(const Tag& t)
    {
        auto it = std::find_if(dirichlet_boundaries_.begin(), dirichlet_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        assert(it != dirichlet_boundaries_.end());
        return *it;
    }

    const MeshCell& Mesh::empty_boundary(const Tag& t) const
    {
        auto it = std::find_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        assert(it != empty_boundaries_.end());
        return *it;
    }

    MeshCell& Mesh::empty_boundary_p(const Tag& t)
    {
        auto it = std::find_if(empty_boundaries_.begin(), empty_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        assert(it != empty_boundaries_.end());
        return *it;
    }

    const MeshCell& Mesh::farfield_boundary(const Tag& t) const
    {
        auto it = std::find_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        assert(it != farfield_boundaries_.end());
        return *it;
    }

    MeshCell& Mesh::farfield_boundary_p(const Tag& t)
    {
        auto it = std::find_if(farfield_boundaries_.begin(), farfield_boundaries_.end(), [&](const MeshCell& tmc){return tmc.tag() == t;});
        assert(it != farfield_boundaries_.end());
        return *it;
    }

    void Mesh::connect_add_bou_to_interior(Mesh& boumesh, boundary_t boutype)
    {
        mcc* container = nullptr;
        if (boutype == "wall")
        {
            container = &wall_boundaries_;
        }
        else if (boutype == "dirichlet")
        {
            container = &dirichlet_boundaries_;
        }
        else if (boutype == "empty")
        {
            container = &empty_boundaries_;
        }
        else if (boutype == "farfield")
        {
            container = &farfield_boundaries_;
        }
        else
        {
            assert(false);
        }
        for (MeshCell mc: boumesh.cell())
        {
            auto pit1 = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return (std::abs(mp.p().r(0) - mc.point(0).p().r(0)) <= ZERO && std::abs(mp.p().r(1) - mc.point(0).p().r(1)) <= ZERO && std::abs(mp.p().r(2) - mc.point(0).p().r(2)) <= ZERO);});
            if (pit1 == point_.end()) continue;

            for (const Tag& pc: pit1->parent_cell())
            {
                auto pit2 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(1).p();});
                if (pit2 != cell(pc).point().end())
                {
                    auto pit3 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(2).p();});
                    if (pit3 != cell(pc).point().end())
                    {
                        //auto pit4 = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(3).p();});

                        mc.set_point_tag(0, pit1->tag());
                        mc.set_point_tag(1, pit2->tag());
                        mc.set_point_tag(2, pit3->tag());

                        if (mc.point().size() > 3)
                        {
                            for (int i=3; i<mc.point().size(); ++i)
                            {
                                auto pit = std::find_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.p() == mc.point(i).p();});
                                mc.set_point_tag(i, pit->tag());
                            }
                        }

                        //assert(cell_p(pc).wall_boundary().size() < 6);
                        cell_p(pc).add_boundary(mc.tag(), boutype);
                        add_element(mc);
                        container->back().set_interior_boundary(pc);

                        MeshCell& mcc = container->back();

                        for (MeshFace& mf: cell_p(pc).face_)
                        {
                            int countt = 0;
                            for (const Tag& mp: mf.mesh_point())
                            {
                                if (mp == pit1->tag() || mp == pit2->tag() || mp == pit3->tag())
                                {
                                    ++countt;
                                }
                            }
                            
                            if (countt >= 3)
                            {
                                mf.add_parent_cell(mcc.tag());
                                mf.add_parent_cell(pc);
                                assert(mcc.face().size() == 1);
                                mcc.face_[0].add_parent_cell(mcc.tag());
                                mcc.face_[0].add_parent_cell(pc);
                                mf.set_as_boundary();
                                mf.set_btype(boutype);
                                mcc.face_[0].set_as_boundary();
                                mcc.face_[0].set_btype(boutype);
                                break;
                            }
                        }

                        break;
                    }
                }
            }
        }
    }

    void Mesh::add_element(const MeshCell& mc)
    {
        if (mc.boutype() == "wall")
        {
            add_boundary(mc, wall_boundaries_, wall_tag_index_map_);
        }
        else if (mc.boutype() == "dirichlet")
        {
            add_boundary(mc, dirichlet_boundaries_, dirichlet_tag_index_map_);
        }
        else if (mc.boutype() == "empty")
        {
            add_boundary(mc, empty_boundaries_, empty_tag_index_map_);
        }
        else if (mc.boutype() == "farfield")
        {
            add_boundary(mc, farfield_boundaries_, farfield_tag_index_map_);
        }
        else if (mc.boutype() == "interior")
        {
            add_interior_cell(mc);
        }
        else
        {
            assert(false);
        }
    }

    /*void Mesh::add_boundary(const MeshCell& mc)
    {
        if (mc.boutype() == boundary_t::wall)
        {
            add_boundary(mc, wall_boundaries_, wall_tag_index_map_);
        }
        else if (mc.boutype() == boundary_t::dirichlet)
        {
            add_boundary(mc, dirichlet_boundaries_, dirichlet_tag_index_map_);
        }
        else if (mc.boutype() == boundary_t::empty)
        {
            add_boundary(mc, empty_boundaries_, empty_tag_index_map_);
        }
        else if (mc.boutype() == boundary_t::farfield)
        {
            add_boundary(mc, farfield_boundaries_, farfield_tag_index_map_);
        }
        else
        {
            assert(false);
        }
    }*/

    void Mesh::add_element_nonsorted(const MeshCell& mc)
    {
        if (mc.boutype() == "wall")
        {
            add_boundary(mc, wall_boundaries_, wall_tag_index_map_);
        }
        else if (mc.boutype() == "dirichlet")
        {
            add_boundary(mc, dirichlet_boundaries_, dirichlet_tag_index_map_);
        }
        else if (mc.boutype() == "empty")
        {
            add_boundary(mc, empty_boundaries_, empty_tag_index_map_);
        }
        else if (mc.boutype() == "farfield")
        {
            add_boundary(mc, farfield_boundaries_, farfield_tag_index_map_);
        }
        else if (mc.boutype() == "interior")
        {
            add_interior_cell_nonsorted(mc);
        }
        else
        {
            assert(false);
        }
    }

    void Mesh::remove_merge_duplicate_points()
    {
        remove_merge_dup(point_);
    }

    void Mesh::add_element_no_check(const MeshCell& mc)
    {
        if (mc.boutype() == "wall")
        {
            add_boundary(mc, wall_boundaries_, wall_tag_index_map_);
        }
        else if (mc.boutype() == "dirichlet")
        {
            add_boundary(mc, dirichlet_boundaries_, dirichlet_tag_index_map_);
        }
        else if (mc.boutype() == "empty")
        {
            add_boundary(mc, empty_boundaries_, empty_tag_index_map_);
        }
        else if (mc.boutype() == "farfield")
        {
            add_boundary(mc, farfield_boundaries_, farfield_tag_index_map_);
        }
        else if (mc.boutype() == "interior")
        {
            add_interior_cell_no_check(mc);
        }
        else
        {
            assert(false);
        }
    }

    /*void Mesh::add_walls(std::vector<MeshCell>& cells)
    {
        assert(wall_boundaries_.size() == wall_tag_index_map_.left.size());
        for (int i=0; i<cells.size(); ++i)
        {
            assert(std::count_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& mc){return mc.tag() == cells[i].tag();}) == 0);
            assert(wall_tag_index_map_.left.count(cells[i].tag()()) == 0);
        }
        //int initial = wall_boundaries_.size();
        for (int i=0; i<cells.size(); ++i)
        {
            if (wall_tag_index_map_.left.count(cells[i].tag()()) != 0) {
                continue;
            }
            assert(wall_tag_index_map_.right.count(wall_boundaries_.size() + i) == 0);
            wall_tag_index_map_.insert(boost::bimap<int, int>::value_type(cells[i].tag()(), wall_boundaries_.size()));
            wall_boundaries_.push_back(cells[i]);
        }
        assert(wall_boundaries_.size() == wall_tag_index_map_.left.size());
    }*/

    /*void Mesh::add_outers(std::vector<MeshCell>& cells)
    {
        assert(outer_boundaries_.size() == outer_tag_index_map_.left.size());
        for (const MeshCell& c: cells)
        {
            assert(std::count_if(outer_boundaries_.begin(), outer_boundaries_.end(), [&](const MeshCell& mc){return mc.tag() == c.tag();}) == 0);
            assert(outer_tag_index_map_.left.count(c.tag()()) == 0);
        }
        //int initial = outer_boundaries_.size();
        for (int i=0; i<cells.size(); ++i)
        {
            if (outer_tag_index_map_.left.count(cells[i].tag()()) != 0) {
                continue;
            }
            assert(outer_tag_index_map_.right.count(outer_boundaries_.size() + i) == 0);
            outer_tag_index_map_.insert(boost::bimap<int, int>::value_type(cells[i].tag()(), outer_boundaries_.size()));
            outer_boundaries_.push_back(cells[i]);
        }
        assert(outer_boundaries_.size() == outer_tag_index_map_.left.size());
    }*/

    void Mesh::add_boundary(const MeshCell& mc, mcc& container, bimap_int& index_map)
    {
        if (std::find_if(container.begin(), container.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != container.end()) {
            return;
        }

        container.push_back(mc);
        index_map.insert(boost::bimap<int, int>::value_type(mc.tag()(), container.size() - 1));
    }

    /*void Mesh::add_wall_boundary(MeshCell&& mc)
    {
        if (std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != wall_boundaries_.end()) {
            return;
        }

        wall_boundaries_.push_back(mc);
        wall_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), wall_boundaries_.size() - 1));
    }

    void Mesh::add_wall_boundary(const MeshCell& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != wall_boundaries_.end()) {
            return;
        }

        wall_boundaries_.push_back(mc);
        wall_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), wall_boundaries_.size() - 1));
    }

    void Mesh::add_outer_boundary(MeshCell&& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(outer_boundaries_.begin(), outer_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != outer_boundaries_.end()) {
            return;
        }
        outer_boundaries_.push_back(mc);
        outer_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), outer_boundaries_.size() - 1));
    }

    void Mesh::add_outer_boundary(const MeshCell& mc)
    {
        // no points are added because interiors, therefore, points are added already.
        if (std::find_if(outer_boundaries_.begin(), outer_boundaries_.end(), [&](const MeshCell& c){return c.tag() == mc.tag();}) != outer_boundaries_.end()) {
            return;
        }

        outer_boundaries_.push_back(mc);
        outer_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), outer_boundaries_.size() - 1));
    }*/

    //void Mesh::add_wall_boundaries(const Mesh& wall_mesh)
    //{
        //for (const MeshCell mc: wall_mesh.cell())
        //{
            //add_wall_boundary(mc);
        //}
    //}

    //void Mesh::add_outer_boundaries(const Mesh& outer_mesh)
    //{
        //for (const MeshCell mc: outer_mesh.cell())
        //{
            //add_outer_boundary(mc);
        //}
    //}

    /*void Mesh::merge_outer_to_interior(const Mesh& outer_mesh)
    {
        for (const MeshCell mc: outer_mesh.cell())
        {
            auto pit = std::find_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.tag() == mc.point(0).tag();});
            if (pit == point_.end()) continue;

            for (const Tag& pc: pit->parent_cell())
            {
                int count = std::count_if(cell(pc).point().begin(), cell(pc).point().end(), [&](const MeshPoint& mp){return mp.tag() == mc.point(1).tag();});
                if (count != 0)
                {
                    cell_p(pc).add_outer_boundary(mc.tag());
                    outer_boundaries_.push_back(mc);
                    outer_boundaries_.back().set_interior_boundary(pc);
                    outer_tag_index_map_.insert(boost::bimap<int, int>::value_type(mc.tag()(), outer_boundaries_.size() - 1));
                    break;
                }
            }
        }
    }*/

    void Mesh::insert_cells(const std::vector<MeshCell>& cells)
    {
        cell_.reserve(cell_.size() + cells.size());
        cell_.insert(cell_.end(), cells.begin(), cells.end());
    }

    void Mesh::sort_cells()
    {
        std::sort(cell_.begin(), cell_.end(), [](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
    }

    void Mesh::sort_points()
    {
        std::sort(point_.begin(), point_.end(), [](const MeshPoint& left, const MeshPoint& right){return left.tag() < right.tag();});
    }

    void Mesh::reserve_interior_only(size_t size)
    {
        cell_.reserve(size);
    }

    void Mesh::reserve_interior(size_t size)
    {
        cell_.reserve(size);
        point_.reserve(size*4); // assuming quad.
    }
    void Mesh::reserve_wall(size_t size)
    {
        wall_boundaries_.reserve(size);
    }
    void Mesh::reserve_farfield(size_t size)
    {
        farfield_boundaries_.reserve(size);
    }
    void Mesh::reserve_empty(size_t size)
    {
        empty_boundaries_.reserve(size);
    }
    void Mesh::reserve_dirichlet(size_t size)
    {
        dirichlet_boundaries_.reserve(size);
    }

    void Mesh::reset_erase_marks()
    {
        for (MeshCell& mc: cell_)
        {
            mc.unmark_to_be_erased();
        }
    }

    void Mesh::mark_to_be_erased(const Tag& ic)
    {
        cell_p(ic).mark_to_be_erased();
    }

    void Mesh::erase_marked_cells()
    {
        for (MeshCell& mc: cell_)
        {
            if (mc.erase())
            {
                prepare_to_remove_cell(mc.tag());
            }
        }
        auto it = std::remove_if(cell_.begin(), cell_.end(), [&](const MeshCell& mc){return mc.erase() == true;});
        cell_.erase(it, cell_.end());
    }

    void Mesh::simple_merge(const Mesh& other_mesh)
    {
        // why interiors are not duplicated but boundaries are?
        assert(tag_ == other_mesh.tag());
        reserve_interior(cell_.size() + other_mesh.cell().size());
        cell_.insert(cell_.end(), other_mesh.cell().begin(), other_mesh.cell().end());
        sort_cells();

        auto aa = [&](const mcc& container)
        {
            for (const MeshCell& mc: container)
            {
                int count = std::count_if(container.begin(), container.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
                if (count != 0) continue;
                add_element(mc);
            }
        };

        aa(other_mesh.wall_boundaries_);
        aa(other_mesh.dirichlet_boundaries_);
        aa(other_mesh.empty_boundaries_);
        aa(other_mesh.farfield_boundaries_);

        /*for (const MeshCell& mc: other_mesh.wall_boundaries_)
        {
            int count = std::count_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            //add_wall_boundary(mc);
            add_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.outer_boundaries_)
        {
            int count = std::count_if(outer_boundaries_.begin(), outer_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            //assert(count == 0);
            if (count != 0) continue;
            add_outer_boundary(mc);
        }*/
    }

    void Mesh::merge_batch(const Mesh& other_mesh, int rank)
    {
        assert(tag_ == other_mesh.tag());
        //reserve_interior(cell_.size() + other_mesh.cell().size());
        cell_.reserve(cell_.size() + other_mesh.cell().size());
        cell_.insert(cell_.end(), other_mesh.cell_.begin(), other_mesh.cell_.end());
        point_.reserve(point_.size() + other_mesh.point().size());
        point_.insert(point_.end(), other_mesh.point_.begin(), other_mesh.point_.end());
        //sort_cells();
        //shrink_cells();
        //shrink_points();
        auto aa = [&](const mcc& container)
        {
            for (const MeshCell& mc: container)
            {
                int count = std::count_if(container.begin(), container.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
                if (count != 0) continue;
                add_element(mc);
            }
        };
        aa(other_mesh.wall_boundaries_);
        aa(other_mesh.dirichlet_boundaries_);
        aa(other_mesh.empty_boundaries_);
        aa(other_mesh.farfield_boundaries_);
    }

    void Mesh::merge(const Mesh& other_mesh)
    {
        assert(tag_ == other_mesh.tag());
        reserve_interior(cell_.size() + other_mesh.cell().size());
        for (const MeshCell& mc: other_mesh.cell())
        {
            //assert(!query(mc.tag()));
            //if (query(mc.tag())) continue;
            //auto it = std::find_if(cell_.begin(), cell_.end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
            //if (it == cell_.end()) {
                //add_interior_cell_sorted(mc);
                add_interior_cell_nonsorted(mc);
            //}
        }
        sort_cells();
        //shrink_cells();
        //shrink_points();
        auto aa = [&](const mcc& container)
        {
            for (const MeshCell& mc: container)
            {
                int count = std::count_if(container.begin(), container.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
                if (count != 0) continue;
                add_element(mc);
            }
        };
        aa(other_mesh.wall_boundaries_);
        aa(other_mesh.dirichlet_boundaries_);
        aa(other_mesh.empty_boundaries_);
        aa(other_mesh.farfield_boundaries_);
    }

    /*void Mesh::make_cell_tags()
    {
        cell_tag_.clear();
        cell_tag_.reserve(cell_.size());
        for (const MeshCell& mc: cell_)
        {
            cell_tag_.push_back(mc.tag()());
        }

        std::sort(cell_tag_.begin(), cell_tag_.end());
    }*/

    /*void Mesh::merge_no_check(const Mesh& other_mesh)
    {
        assert(tag_ == other_mesh.tag());
        reserve_interior(cell_.size() + other_mesh.cell().size());
        for (const MeshCell& mc: other_mesh.cell())
        {
            add_interior_cell_no_check(mc);
        }
        for (const MeshCell& mc: other_mesh.wall_boundaries_)
        {
            auto it = std::find_if(wall_boundaries_.begin(), wall_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            if (it != wall_boundaries_.end()) continue;
            add_wall_boundary(mc);
        }
        for (const MeshCell& mc: other_mesh.outer_boundaries_)
        {
            auto it = std::find_if(outer_boundaries_.begin(), outer_boundaries_.end(), [&](const MeshCell& wb){return wb.tag() == mc.tag();});
            if (it != outer_boundaries_.end()) continue;
            add_outer_boundary(mc);
        }
        remove_merge_dup_pts();
        remove_dup_cells();
        set_point_tag_index_map();
        set_cell_tag_index_map();
        shrink_points();
        shrink_cells();
    }*/

    void Mesh::remove_merge_duplicate_cells()
    {
        remove_merge_dup(cell_);
    }

    void Mesh::remove_duplicate_cells()
    {
        remove_dup(cell_);
    }

    void Mesh::remove_duplicate_points()
    {
        remove_dup(point_);
    }

    //void Mesh::receptor_to_field(std::deque<Mesh>& mesh)
    //{
        //for (MeshCell& mc: cell_)
        //{
            //if (mc.oga_cell_type() != OGA_cell_type_t::mandat_receptor) continue;
            ////assert(mc.donor_mesh().isvalid());
            ////assert(mc.donor_cell().isvalid());
            //assert(mc.donor().first.isvalid());
            //assert(mc.donor().second.isvalid());

            //for (Mesh& m: mesh)
            //{
                ////if (m.tag() != mc.donor_mesh()) continue;
                //if (m.tag() != mc.donor().first) continue;

                ////m.cell_p(mc.donor_cell()).set_oga_cell_type(OGA_cell_type_t::field);
                //m.cell_p(mc.donor().second).set_oga_cell_type(OGA_cell_type_t::field);
                //break;
            //}
        //}
    //}

    //void Mesh::fringe_to_field()
    //{
        //for (MeshCell& mc: cell_)
        //{
            ////if (mc.is_ghost()) continue;

            ////if (mc.polygon().edge().size() < mc.nei().size())
            ////{
            ////std::cout << mc.polygon().edge().size() << std::endl;
            ////std::cout << mc.nei().size() << std::endl;
            ////}
            ////assert(mc.polygon().edge().size() >= mc.nei().size());
            ////assert(!mc.pnei().empty());

            //if (mc.oga_cell_type() != OGA_cell_type_t::receptor) continue;
            //if (mc.oga_cell_type() == OGA_cell_type_t::mandat_receptor) continue;

            //int count = 0;
            //for (const Tag& n: mc.pnei())
            //{
                //auto it = cell_tag_index_map.left.find(n());
                //if (it == cell_tag_index_map.left.end())
                //{
                    //std::cout << mc.pnei().size() << std::endl;
                    //std::cout << n() << std::endl;
                //}
                //assert(it != cell_tag_index_map.left.end());

                ////if (cell(n).is_ghost()) continue;

                //if (cell(n).oga_cell_type() == OGA_cell_type_t::field) ++count;
            //}

            //if (count != 0 && count == mc.pnei().size())
                //mc.set_oga_cell_type(OGA_cell_type_t::field);
        //}
    //}

    const MeshPoint* Mesh::query_point(const Tag& ic) const
    {
        assert(ic.isvalid());

        auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ic, Tag()));
        if (it != point_.end() && it->tag() == ic)
        {
            return &*it;
        }

        return nullptr;
    }

    const MeshCell* Mesh::query_sorted(const Tag& ic) const
    {
        assert(ic.isvalid());

        MeshCell tmc;
        tmc.set_tag(ic);

        auto it = std::lower_bound(cell_.begin(), cell_.end(), tmc, [](const MeshCell& left, const MeshCell& right){return left.tag()() < right.tag()();});

        if (it != cell_.end())
        {
            if (it->tag() == ic) {
                return &(*it);
            }
        }

        return nullptr;
    }

    const MeshCell* Mesh::query(const Tag& ic) const
    {
        assert(ic.isvalid());

        auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ic, Tag()));
        if (it != cell_.end() && it->tag() == ic)
        {
            return &*it;
        }

        return nullptr;
    }

    const MeshCell* Mesh::query_wall(const Tag& ic) const
    {
        assert(ic.isvalid());

        int count = wall_tag_index_map_.left.count(ic());
        //auto it = cell_tag_index_map.left.find(ic());

        if (count != 0) {
        //if (it != cell_tag_index_map.left.end()) {
            return &wall_boundary(ic);
        }

        //if (count != 0) return true;

        return nullptr;
    }

    const MeshCell* Mesh::query_dirichlet(const Tag& ic) const
    {
        assert(ic.isvalid());

        int count = dirichlet_tag_index_map_.left.count(ic());
        //auto it = cell_tag_index_map.left.find(ic());

        if (count != 0) {
        //if (it != cell_tag_index_map.left.end()) {
            return &dirichlet_boundary(ic);
        }

        //if (count != 0) return true;

        return nullptr;
    }

    const AABB& Mesh::hole_aabb() const
    {
        return hole_aabb_;
    }

    void Mesh::remove_cell_from_cellhood(const Tag& ic)
    {
        assert(ic.isvalid());
        assert(query(ic));
        MeshCell& mc = cell_p(ic);

        for (const Tag& inei: mc.pnei())
        {
            assert(inei.isvalid());
            assert(query(inei));
            cell_p(inei).remove_neighbor(ic);
            //std::find(cell(inei).pnei().begin(), cell(inei).pnei().end(), )cell(inei).pnei().find();
        }

        mc.remove_all_neighbors();
    }

    void Mesh::deparent_cell_from_vertices(const Tag& ic)
    {
        for (const MeshPoint& mp: cell_p(ic).point())
        {
            point_p(mp.tag()).remove_parent_cell(ic);
        }
    }

    unsigned int Mesh::point_index(const Tag& t) const
    {
        auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(t, Tag()));
        assert(it != point_.end());
        assert(it->tag() == t);
        return std::distance(point_.begin(), it);
    }

    void Mesh::rotate(double angle, int axis, const vec3<double>& rot_point)
    {
        // rotate points.
        for (MeshPoint& point: point_)
        {
            point.rotate_point(angle, axis, rot_point);
        }
        // rotate cell points.
        for (MeshCell& cell: cell_)
        {
            cell.rotate_points(angle, axis, rot_point);
        }
        auto aa = [&](mcc& container)
        {
            for (MeshCell& cell: container)
            {
                cell.rotate_points(angle, axis, rot_point);
            }
        };
        // rotate cell points.
        aa(wall_boundaries_);
        aa(dirichlet_boundaries_);
        aa(empty_boundaries_);
        aa(farfield_boundaries_);
    }

    void Mesh::move(const vec3<double>& v)
    {
        // this is a temporary function to mimic solver.
        // aim is to displace certain meshblocks.

        {
            // move points.
            for (MeshPoint& point: point_)
            {
                point.move_point(v);
            }
            // move cell points.
            for (MeshCell& cell: cell_)
            {
                cell.move_points(v);
            }

            auto aa = [&](mcc& container)
            {
                for (MeshCell& cell: container)
                {
                    cell.move_points(v);
                }
            };

            aa(wall_boundaries_);
            aa(dirichlet_boundaries_);
            aa(empty_boundaries_);
            aa(farfield_boundaries_);
        }
    }

    void Mesh::set_cell(const Tag& t, const MeshCell& c)
    {
        cell_p(t) = c;
    }

    const std::vector<MeshPoint>& Mesh::point() const
    {
        return point_;
    }

    std::vector<Point> Mesh::rawpoint() const
    {
        std::vector<Point> rawpoints;
        rawpoints.reserve(point_.size());
        for (const auto& p: point_)
        {
            rawpoints.push_back(p.p());
        }

        return rawpoints;
    }

    const Tag& Mesh::tag() const
    {
        return tag_;
    }

    const Tag& Mesh::parent_mesh() const
    {
        return parent_mesh_;
    }

    void Mesh::set_tag(const Tag& t)
    {
        tag_ = t;
    }
    void Mesh::set_parent_mesh(const Tag& t)
    {
        parent_mesh_ = t;
    }

    const std::vector<MeshCell>& Mesh::cell() const
    {
        return cell_;
    }

    /*std::vector<MeshCell>& Mesh::cell()
      {
      return cell_;
      }*/

    MeshPoint& Mesh::point_p(const Tag& ptag)
    {
        assert(ptag.isvalid());

        auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ptag, Tag()));
        assert(it != point_.end());
        assert(it->tag() == ptag);

        return *it;
    }

    const MeshPoint& Mesh::point(const Tag& ptag) const
    {
        assert(ptag.isvalid());

        auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ptag, Tag()));
        assert(it != point_.end());
        assert(it->tag() == ptag);

        return *it;
    }

    const MeshCell* Mesh::cell_ptr(const Tag& ctag) const
    {
        assert(ctag.isvalid());

        auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ctag, Tag()));
        if (it == cell_.end())
        {
            return nullptr;
        }
        else if (it->tag() != ctag)
        {
            return nullptr;
        }

        return &(*it);
    }

    MeshFace& Mesh::face_p(const Tag& ctag)
    {
        assert(ctag.isvalid());

        auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& mf){return mf.tag() == ctag;});
        assert(it != face_.end());
        if (it->tag() != ctag)
        {
            std::cout << "it->tag(): " << it->tag()() << std::endl;
            std::cout << "ctag: " << ctag() << std::endl;
        }
        assert(it->tag() == ctag);

        return *it;
    }
    
    const MeshFace& Mesh::face(const Tag& ctag) const
    {
        assert(ctag.isvalid());

        auto it = std::find_if(face_.begin(), face_.end(), [&](const MeshFace& mf){return mf.tag() == ctag;});
        assert(it != face_.end());
        if (it->tag() != ctag)
        {
            std::cout << "it->tag(): " << it->tag()() << std::endl;
            std::cout << "ctag: " << ctag() << std::endl;
        }
        assert(it->tag() == ctag);

        return *it;
    }

    const MeshCell& Mesh::cell(const Tag& ctag) const
    {
        assert(ctag.isvalid());

        auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ctag, Tag()));
        assert(it != cell_.end());
        if (it->tag() != ctag)
        {
            std::cout << "it->tag(): " << it->tag()() << std::endl;
            std::cout << "ctag: " << ctag() << std::endl;
        }
        assert(it->tag() == ctag);

        return *it;
    }

    MeshCell& Mesh::cell_p(const Tag& ctag)
    {
        assert(ctag.isvalid());

        auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ctag, Tag()));
        assert(it != cell_.end());
        assert(it->tag() == ctag);

        return *it;
    }

    void Mesh::add_cell_only(MeshCell c, size_t size)
    {
        auto add = [&] ()
        {
            /*c.set_parent_cell_of_vertices();
            assert(tag().isvalid());
            for (const MeshPoint& p: c.point())
            {
                assert(c.tag().isvalid());
                point_p(p.tag()).add_parent_cell(c.tag());
            }
            // add to container.
            cell_.push_back(std::move(c));
            for (const MeshPoint& p: cell_.back().point())
            {
                assert(p.parent_cell().size() > 0);
            }
            cell_tag_index_map.insert(boost::bimap<int, int>::value_type(c.tag()(), cell_.size() - 1));*/

            bool dummy;
            //add_cell_sorted(c, dummy);
            add_cell_nonsorted(c);
        };

        // reserve memory for cell_ if size is provided.
        if (size != 0)
        {
            if (cell_.capacity() != size)
            {
                cell_.reserve(size);
            }
        }

        assert(c.tag().isvalid());
        if (c.tag().isvalid())
        {
            //auto it = std::find_if(cell_.begin(), cell_.end(), [&c](const MeshCell& _c){return _c.tag() == c.tag();});
            //assert(it == cell_.end());
            //if (it == cell_.end())
            {
                add();
            }
        }
        else
        {
            add();
        }
    }

    void Mesh::shrink_cells()
    {
        mcc(cell_).swap(cell_);
    }

    void Mesh::shrink_points()
    {
        mpc(point_).swap(point_);
    }

    /*void Mesh::add_interior_cell(MeshCell&& c)
    {
        c.set_parent_mesh(tag());
        c.remove_parent_cells_of_vertices();
        cell_.push_back(c);
        cell_tag_index_map.insert(boost::bimap<int, int>::value_type(cell_.back().tag()(), cell_.size() - 1));
        cell_.back().set_parent_cell_of_vertices();

        for (const MeshPoint& p: cell_.back().point())
        {
            auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
            if (it == point_.end())
            { 
                point_.push_back(p);
                //auto itt = point_tag_index_map.left.find(p.tag()());
                //assert(itt == point_tag_index_map.left.end());
                point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
            }
            else
            {
                it->add_parent_cell(cell_.back().tag());
            }
        }
    }*/

    void Mesh::remove_dup_cells_and_points()
    {
        ////remove_merge_dup_pts();
        remove_merge_duplicate_cells(); // because there may be dup cells after add and push sp in remap. dups are due to being right on border of two or more bins.
        remove_merge_duplicate_points();

        ////set_point_tag_index_map();
        ////set_cell_tag_index_map();
        
        shrink_points();
        shrink_cells();
    }

    void Mesh::add_interior_nonsorted_addpoint(MeshCell mc)
    {
        mc.set_parent_mesh(tag());
        for (const MeshPoint& p: mc.point())
        {
            point_.push_back(p);
        }
        cell_.push_back(std::move(mc));
    }

    void Mesh::add_interiors_nonsorted_nopoint(std::vector<MeshCell>& cells)
    {
        for (MeshCell& mc: cells)
        {
            mc.set_parent_mesh(tag());
            //mc.remove_parent_cells_of_vertices();
            //mc.set_parent_cell_of_vertices();
        }

        cell_.insert(cell_.end(), cells.begin(), cells.end());
    }

    void Mesh::add_interior_cell_no_check(const MeshCell& c)
    {
        assert(tag().isvalid());
        assert(c.tag().isvalid());
        //assert(!query(c.tag()));

        cell_.push_back(c);
        cell_.back().set_parent_mesh(tag());
        cell_.back().remove_parent_cells_of_vertices();

        //assert(query(cell_.back().tag()));

        cell_.back().set_parent_cell_of_vertices();

        assert(!cell_.back().point().empty());

        point_.insert(point_.end(), cell_.back().point().begin(), cell_.back().point().end());

        for (const MeshPoint& p: cell_.back().point())
        {
            assert(!p.parent_cell().empty());
        }

        for (const MeshPoint& _p: point_)
        {
            for (const Tag& _t: _p.parent_cell())
            {
                assert(_t.isvalid());
            }
        }
    }

    /*void Mesh::set_cell_tag_index_map()
    {
        cell_tag_index_map.clear();
        for (const MeshCell& mc: cell_)
        {
            cell_tag_index_map.insert(boost::bimap<int, int>::value_type(mc.tag()(), cell_tag_index_map.size()));
        }
    }

    void Mesh::set_point_tag_index_map()
    {
        point_tag_index_map.clear();
        for (const MeshPoint& mp: point_)
        {
            point_tag_index_map.insert(boost::bimap<int, int>::value_type(mp.tag()(), point_tag_index_map.size()));
        }
    }*/

    /*void Mesh::add_interior_cell_find(const MeshCell& c)
    {
        assert(tag().isvalid());
        assert(c.tag().isvalid());
        assert(!query(c.tag()));

        cell_.push_back(c);
        cell_.back().set_parent_mesh(tag());
        cell_.back().remove_parent_cells_of_vertices();
        cell_tag_index_map.insert(boost::bimap<int, int>::value_type(c.tag()(), cell_.size() - 1));

        assert(query(cell_.back().tag()));

        cell_.back().set_parent_cell_of_vertices();

        assert(!cell_.back().point().empty());
        
        for (const MeshPoint& p: cell_.back().point())
        {
            //assert(p.parent_mesh() == c.parent_mesh());

            auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
            if (it == point_.end())
            { 
                assert(cell_.back().tag() == c.tag());
                assert(p.tag().isvalid());
                point_.push_back(p);
                //auto itt = point_tag_index_map.left.find(p.tag()()); 
                //assert(itt == point_tag_index_map.left.end());
                point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                //for (const Tag& _t: point_.back().parent_cell())
                //{
                    //assert(_t.isvalid());
                //}

            }
            else
            {
                it->add_parent_cell(cell_.back().tag());
            }
        }

        for (const MeshPoint& p: cell_.back().point())
        {
            assert(!p.parent_cell().empty());
        }

        for (const MeshPoint& _p: point_)
        {
            for (const Tag& _t: _p.parent_cell())
            {
                assert(_t.isvalid());
            }
        }
    }*/

    std::vector<MeshCell>::iterator Mesh::add_cell_sorted(const MeshCell& c, bool& exist)
    {
        auto it = std::lower_bound(cell_.begin(), cell_.end(), c, [](const MeshCell& left, const MeshCell& right){return left.tag()() < right.tag()();});

        if (it == cell_.end())
        { 
            cell_.push_back(c);
            it = std::prev(cell_.end());
        }
        else if (it->tag() != c.tag())
        {
            //assert(it->tag() != c.tag());
            it = cell_.insert(it, c);
        }
        else
        {
            exist = true;
            return it;
        }

        it->set_parent_mesh(tag_);
        it->remove_parent_cells_of_vertices();
        it->set_parent_cell_of_vertices();
        
        assert(!it->point().empty());

        return it;
    }

    void Mesh::remove_parent_cells_of_all_points()
    {
        for (MeshPoint& mp: point_)
        {
            mp.remove_parent_cells();
        }
    }

    void Mesh::remove_parent_cells_of_vertices_of_all_cells()
    {
        for (MeshCell& mc: cell_)
        {
            mc.remove_parent_cells_of_vertices();
        }
    }

    void Mesh::set_parent_cell_of_vertices_of_all_cells()
    {
        for (MeshCell& mc: cell_)
        {
            mc.set_parent_cell_of_vertices();
        }

        for (const MeshPoint& mp: point())
        {
            for (const Tag& pc: mp.parent_cell())
            {
                assert(query(pc) != nullptr);
            }
        }
    }

    void Mesh::update_points_from_cell_vertices()
    {
        // do this after points are sorted and mesh cells vertices are ready.
        for (const MeshCell& mc: cell_)
        {
            for (const MeshPoint& mp: mc.point())
            {
                auto it = std::lower_bound(point_.begin(), point_.end(), mp, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
                assert(it != point_.end());
                assert(it->tag() == mp.tag());
                it->add_parent_cell(mc.tag());

                for (const Tag& pc: it->parent_cell())
                {
                    assert(query(pc) != nullptr);
                }
            }
        }
    }

    std::vector<MeshCell>::iterator Mesh::add_cell_nonsorted(const MeshCell& c)
    {
        cell_.push_back(std::move(c));
        auto it = std::prev(cell_.end());
        
        //it->set_parent_mesh(tag_);
        it->remove_parent_cells_of_vertices();
        it->set_parent_cell_of_vertices();
        
        assert(!it->point().empty());

        return it;
    }

    void Mesh::add_points(const std::vector<MeshCell>& cells)
    {
        std::vector<const MeshPoint*> addresses;
        for (const MeshCell& mc: cells)
        {
            for (const MeshPoint& p: mc.point())
            {
                auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});

                if (it == point_.end())
                { 
                    assert(p.tag().isvalid());
                    point_.push_back(p);
                }
                else if (it->tag() == p.tag())
                {
                    it->add_parent_cell(mc.tag());
                }
                else
                {
                    point_.insert(it, p);
                }
            }
        }

        std::sort(point_.begin(), point_.end());
    }

    void Mesh::add_interior_cell_sorted(const MeshCell& c)
    {
        assert(tag().isvalid());
        assert(c.tag().isvalid());
        //assert(!query(c.tag()));

        bool exist = false;
        auto mc = add_cell_sorted(c, exist);
        if (exist) {
            return;
        }
        
        for (const MeshPoint& p: mc->point())
        {
            //assert(p.parent_mesh() == c.parent_mesh());

            //auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
            // if point_ is not empty, must be already sorted.
            auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
            if (it == point_.end())
            { 
                assert(mc->tag() == c.tag());
                assert(p.tag().isvalid());
                point_.push_back(p);
                //auto itt = point_tag_index_map.left.find(p.tag()()); 
                //assert(itt == point_tag_index_map.left.end());
                //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                //for (const Tag& _t: point_.back().parent_cell())
                //{
                    //assert(_t.isvalid());
                //}

            }
            else if (it->tag() == p.tag())
            {
                it->add_parent_cell(mc->tag());
            }
            else
            {
                point_.insert(it, p);
            }
        }

        for (const MeshPoint& p: cell_.back().point())
        {
            assert(!p.parent_cell().empty());
        }

        for (const MeshPoint& _p: point_)
        {
            for (const Tag& _t: _p.parent_cell())
            {
                assert(_t.isvalid());
            }
        }
    }

    //void Mesh::add_meshpoint(const MeshCell& mc)
    void Mesh::add_meshpoint(const MeshCell& mc, int rank)
    {
        for (const MeshPoint& p: mc.point())
        {
            if (rank == 6) {std::cout << "before lower" << std::endl;}
            auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
            if (rank == 6) {std::cout << "after lower" << std::endl;}
            if (it == point_.end())
            { 
                if (rank == 6) {std::cout << "admmeshpoint0" << std::endl;}
                assert(p.tag().isvalid());
                point_.push_back(p);

            }
            else if (it->tag() == p.tag())
            {
                if (rank == 6) {std::cout << "admmeshpoint1" << std::endl;}
                it->add_parent_cell(mc.tag());
            }
            else
            {
                if (rank == 6) {std::cout << "admmeshpoint2" << std::endl;}
                if (rank == 6) {std::cout << "it tag = " << it->tag()() << std::endl;}
                if (rank == 6) {std::cout << "p tag = " << p.tag()() << std::endl;}
                if (rank == 6) {std::cout << "point size = " << point_.size() << std::endl;}
                point_.insert(it, p);
                if (rank == 6) {std::cout << "admmeshpoint3" << std::endl;}
            }
            if (rank == 6) {std::cout << "admmeshpoint4" << std::endl;}
        }
        if (rank == 6) {std::cout << "admmeshpoint5" << std::endl;}
    }

    void Mesh::add_interior_cell_nonsorted(const MeshCell& c)
    {
        assert(tag().isvalid());
        assert(c.tag().isvalid());
        //assert(!query(c.tag()));

        auto mc = add_cell_nonsorted(c);
        
        for (const MeshPoint& p: mc->point())
        {
            //assert(p.parent_mesh() == c.parent_mesh());

            //auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
            // if point_ is not empty, must be already sorted.
            auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
            if (it == point_.end())
            { 
                assert(mc->tag() == c.tag());
                assert(p.tag().isvalid());
                point_.push_back(p);
                //auto itt = point_tag_index_map.left.find(p.tag()()); 
                //assert(itt == point_tag_index_map.left.end());
                //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                //for (const Tag& _t: point_.back().parent_cell())
                //{
                    //assert(_t.isvalid());
                //}

            }
            else if (it->tag() == p.tag())
            {
                it->add_parent_cell(mc->tag());
            }
            else
            {
                point_.insert(it, p);
            }
        }

        for (const MeshPoint& p: cell_.back().point())
        {
            assert(!p.parent_cell().empty());
        }

        for (const MeshPoint& _p: point_)
        {
            for (const Tag& _t: _p.parent_cell())
            {
                assert(_t.isvalid());
            }
        }
    }

    void Mesh::add_interior_cell(const MeshCell& c)
    {
        assert(tag().isvalid());
        assert(c.tag().isvalid());
        assert(!query(c.tag()));

        //cell_.push_back(c);
        //cell_.back().set_parent_mesh(tag());
        //cell_.back().remove_parent_cells_of_vertices();
        //cell_tag_index_map.insert(boost::bimap<int, int>::value_type(c.tag()(), cell_.size() - 1));

        //assert(query(cell_.back().tag()));

        //cell_.back().set_parent_cell_of_vertices();

        //assert(!cell_.back().point().empty());

        bool dummy = false;
        auto mc = add_cell_sorted(c, dummy);
        
        for (const MeshPoint& p: mc->point())
        {
            //assert(p.parent_mesh() == c.parent_mesh());

            //auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& mp){return mp.tag() == p.tag();});
            // if point_ is not empty, must be already sorted.
            auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
            if (it == point_.end())
            { 
                assert(mc->tag() == c.tag());
                assert(p.tag().isvalid());
                point_.push_back(p);
                //auto itt = point_tag_index_map.left.find(p.tag()()); 
                //assert(itt == point_tag_index_map.left.end());
                //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
                //for (const Tag& _t: point_.back().parent_cell())
                //{
                    //assert(_t.isvalid());
                //}

            }
            else if (it->tag() == p.tag())
            {
                it->add_parent_cell(mc->tag());
            }
            else
            {
                point_.insert(it, p);
            }
        }

        for (const MeshPoint& p: cell_.back().point())
        {
            assert(!p.parent_cell().empty());
        }

        for (const MeshPoint& _p: point_)
        {
            for (const Tag& _t: _p.parent_cell())
            {
                assert(_t.isvalid());
            }
        }
    }

    void Mesh::add_point(MeshPoint p, const Tag& t, size_t size)
    {
        // reserve memory for cell_ if size is provided.
        if (size != 0)
        {
            if (point_.capacity() != size)
            {
                point_.reserve(size);
            }
        }

        //auto it = std::find_if(point_.begin(), point_.end(), [&p](const MeshPoint& _p){return _p.tag() == p.tag();});
        //auto it = std::lower_bound(point_.begin(), point_.end(), p, [](const MeshPoint& left, const MeshPoint& right){return left.tag()() < right.tag()();});
        //assert(it == point_.end());
        //if (it == point_.end())
        {
            // generate a new cell tag.
            //p.set_tag(pointtag);
            p.set_tag(t);
            // set parent mesh.
            p.set_parent_mesh(tag_);
            // add to container.
            point_.push_back(std::move(p));
            // insert to bimap.
            //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
        }
        /*else
        {
            p.set_tag(t);
            p.set_parent_mesh(tag_);

            assert(it->tag() != p.tag());
            point_.insert(it, p);
            //point_tag_index_map.insert(boost::bimap<int, int>::value_type(p.tag()(), point_.size() - 1));
        }*/
    }

    bool Mesh::do_point_exist(const Tag& t) const
    {
        auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(t, Tag()));
        if (it != point_.end() && it->tag() == t)
        {
            return true;
        }

        return false;
    }

    void Mesh::remove_cell_boundary(const Tag& ic)
    {
        assert(dirichlet_boundaries_.size() == dirichlet_tag_index_map_.size());
        const MeshCell& mc = cell(ic);

        auto aa = [&](const MeshCell::bouarray& mccontainer, mcc& container, bimap_int& index_map)
        {
            for (const Tag& t: mccontainer)
            {
                auto it = index_map.left.find(t());
                assert(it != index_map.left.end());
                container.erase(container.begin() + it->second);

                int thres = it->second;
                index_map.left.erase(t());
                for (int j=thres+1; j<=container.size(); ++j)
                {
                    auto itt = index_map.right.find(j);
                    if (itt == index_map.right.end())
                    {
                        std::cout << "size map: " << index_map.size() << std::endl;
                        std::cout << "size vec: " << container.size() << std::endl;
                        std::cout << "j: " << j << std::endl;
                        std::cout << "thres: " << thres << std::endl;
                        for (auto yy=index_map.right.begin(); yy != index_map.right.end(); ++yy)
                        {
                            std::cout << "yy: " << yy->first << std::endl;
                        }
                    }
                    assert(itt != index_map.right.end());
                    int rep = itt->first - 1;
                    bool successful_replace = index_map.right.replace_key(itt, rep);
                    assert(successful_replace);
                    if (container[rep].tag()() != itt->second)
                    {
                        std::cout << "wb size: " << container.size() << std::endl;
                        std::cout << "map size: " << index_map.size() << std::endl;
                        std::cout << "rep: " << container[rep].tag()() << std::endl;
                        std::cout << "second: " << itt->second << std::endl;
                        std::cout << "j: " << j << std::endl;
                        std::cout << "thres: " << thres << std::endl;
                    }
                    assert(container[rep].tag()() == itt->second);
                }
                assert(index_map.left.count(t()) == 0);
            }
        };

        aa(mc.wall_boundary(), wall_boundaries_, wall_tag_index_map_);
        aa(mc.dirichlet_boundary(), dirichlet_boundaries_, dirichlet_tag_index_map_);
        aa(mc.empty_boundary(), empty_boundaries_, empty_tag_index_map_);
        aa(mc.farfield_boundary(), farfield_boundaries_, farfield_tag_index_map_);
    }

    const MeshCell& Mesh::boundary(boundary_t btype, const Tag& t) const
    {
        if (btype == "wall")
        {
            return wall_boundary(t);
        }
        else if (btype == "dirichlet")
        {
            return dirichlet_boundary(t);
        }
        else if (btype == "farfield")
        {
            return farfield_boundary(t);
        }
        else if (btype == "empty")
        {
            return empty_boundary(t);
        }
        else
        {
            assert(false);
        }
    }

    MeshCell& Mesh::boundary(boundary_t btype, const Tag& t)
    {
        if (btype == "wall")
        {
            return wall_boundary_p(t);
        }
        else if (btype == "dirichlet")
        {
            return dirichlet_boundary_p(t);
        }
        else if (btype == "farfield")
        {
            return farfield_boundary_p(t);
        }
        else if (btype == "empty")
        {
            return empty_boundary_p(t);
        }
        else
        {
            assert(false);
        }
    }

    void Mesh::prepare_to_remove_cell(const Tag& ic)
    {
        assert(ic.isvalid());

        deparent_cell_from_vertices(ic);
        cell_p(ic).deparent_self_from_faces();
        
        // remove self from neighbors' faces as well.
        for (const Tag& inei: cell_p(ic).pnei())
        {
            cell_p(inei).deparent_cell_from_faces(ic);
        }

        remove_cell_from_cellhood(ic); // needed to update neighbors.
        remove_cell_boundary(ic);
    }

    void Mesh::remove_cell(const Tag& ic)
    {
        prepare_to_remove_cell(ic);

        // erase cell.
        //auto it = cell_tag_index_map.left.find(ic());
        //assert(it != cell_tag_index_map.left.end());
        //cell_.erase(cell_.begin() + it->second);
        
        auto it = std::lower_bound(cell_.begin(), cell_.end(), MeshCell(ic, Tag()));
        assert(it != cell_.end() && it->tag() == ic);
        cell_.erase(it);

        // remove from cell map.
        // example:
        // a b c d e f g h
        // 0 1 2 3 4 5 6 7
        // 0 1 2 x 4 5 6 7
        // thres+1 = 4, cell_.size()=7
        /*int thres = it->second;
        cell_tag_index_map.left.erase(ic());
        for (int j=thres+1; j<=cell_.size(); ++j)
        {
            auto itt = cell_tag_index_map.right.find(j);
            assert(itt != cell_tag_index_map.right.end());
            int rep = itt->first - 1;
            bool successful_replace = cell_tag_index_map.right.replace_key(itt, rep);
            assert(successful_replace);
            assert(cell_[rep].tag()() == itt->second);
            //assert(query(itt->second));
        }
        assert(cell_tag_index_map.left.count(ic()) == 0);*/

    }

    void Mesh::remove_isolated_points()
    {
        for (MeshPoint& mp: point_)
        {
            if (mp.parent_cell().empty())
            {
                mp.mark_to_be_erased();
            }
        }

        auto it = std::remove_if(point_.begin(), point_.end(), [&](const MeshPoint& mp){return mp.erase() == true;});
        point_.erase(it, point_.end());
    }

    void Mesh::remove_point(Tag ip)
    {
        assert(ip.isvalid());

        auto it = std::lower_bound(point_.begin(), point_.end(), MeshPoint(ip, Tag()));
        assert(it != point_.end() && it->tag() == ip);

        point_.erase(it);

        /*// remove from map.
        int thres = it->second;
        point_tag_index_map.left.erase(ip());
        for (int j=thres+1; j<=point_.size(); ++j)
        {
            auto itt = point_tag_index_map.right.find(j);
            assert(itt != point_tag_index_map.right.end()); 
            int rep = itt->first - 1;
            bool successful_replace = point_tag_index_map.right.replace_key(itt, rep);
            assert(successful_replace);
        }*/
    }

    size_t Mesh::calc_nhole_cell() const
    {
        size_t nhole_cell = 0;
        for (const MeshCell& mc: cell_)
        {
            if (mc.oga_cell_type() == OGA_cell_type_t::hole)
                ++nhole_cell;
        }

        return nhole_cell;
    }

    Mesh::Mesh(const Tag& tag, const Tag& parent_mesh): tag_(tag), parent_mesh_(parent_mesh)
    {
        //if (file_name != "")
        //{
            //if (rank != 0)
            //{
                //file_name.append("_");
                //file_name.append(std::to_string(rank));
            //}
            //file_name.append(".msh");
            ////std::cout << "file name: " << file_name << std::endl;
            //read_mesh_GMSH(file_name);
        //}

        //hole_map_.set_nstripe(1, 1);

        //determine_hole_aabb();
    }
}
