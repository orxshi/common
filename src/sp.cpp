#include "sp.h"

namespace Common
{
    const std::deque<ADT>& SpatialPartition::cell_adt() const
    {
        return cell_adt_;
    }

    const std::deque<AABB>& SpatialPartition::mesh_aabb() const
    {
        return mesh_aabb_;
    }

    const AABB& SpatialPartition::mesh_aabb(int i) const
    {
        return mesh_aabb_[i];
    }

    const ADT& SpatialPartition::cell_adt(int i) const
    {
        return cell_adt_[i];
    }

    void SpatialPartition::make_mesh_aabb()
    {
        for (const Mesh& m: mesh_)
        {
            mesh_aabb_.push_back(AABB(m.rawpoint()));
        }
    }

    void SpatialPartition::make_cell_adt()
    {
        for (const Mesh& m: mesh_)
        {
            std::vector<ADTPoint> adt_pts;
            adt_pts.reserve(m.cell().size());
            for (const MeshCell& mc: m.cell())
            {
                std::vector<Point> pts;
                pts.reserve(mc.point().size());
                for (auto it=mc.point().begin(); it != mc.point().end(); ++it) {
                    pts.push_back(it->p());
                }
                adt_pts.push_back(ADTPoint(pts.begin(), pts.end(), mc.tag()()));
            }

            if (adt_pts.empty())
            {
                cell_adt_.push_back(ADT());
            }
            else
            {
                cell_adt_.push_back(ADT(adt_pts));
            }
        }
    }

    void SpatialPartition::extend_aabb(const AABB& aabb)
    {
        aabb_.extend(aabb);
    }

    void SpatialPartition::profstart(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->start(fun);
    }

    void SpatialPartition::profstop(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->stop(fun);
    }

    void SpatialPartition::set_profiler(std::shared_ptr<Profiler> profiler)
    {
        profiler_ = profiler;
    }

    void SpatialPartition::merge_meshes_no_check(std::deque<Mesh>& mesh) const
    {
        for (const Mesh& m: mesh_)
        {
            auto mit = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& tmpm){return tmpm.tag() == m.tag();});
            
            if (mit != mesh.end())
            {
                //mit->reserve_interior(m.cell().size()); // wrong. should be current size + incoming size.
                mit->reserve_interior(mit->cell().size() + m.cell().size());
                for (const MeshCell& mc: m.cell())
                {
                    mit->add_element_no_check(mc);
                }
                auto aa = [&](const mcc& container)
                {
                    for (const MeshCell& mc: container)
                    {
                        int count = std::count_if(container.begin(), container.end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                        if (count == 0)
                        {
                            mit->add_element(mc);
                        }
                    }
                };
                aa(m.wall_boundaries());
                aa(m.dirichlet_boundaries());
                aa(m.empty_boundaries());
                aa(m.farfield_boundaries());
            }
            else
            {
                mesh.push_front(m);
            }
        }
    }

    void SpatialPartition::merge_meshes(std::deque<Mesh>& mesh) const
    {
        assert(!mesh_.empty());
        for (const Mesh& m: mesh_)
        {
            auto mit = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& tmpm){return tmpm.tag() == m.tag();});
            
            if (mit != mesh.end())
            {
                for (const MeshCell& mc: m.cell())
                {
                    int count = std::count_if(mit->cell().begin(), mit->cell().end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                    if (count == 0)
                    {
                        mit->add_element_nonsorted(mc);
                    }
                }
                auto aa = [&](const mcc& container)
                {
                    for (const MeshCell& mc: container)
                    {
                        int count = std::count_if(container.begin(), container.end(), [&](const MeshCell& tmc){return tmc.tag() == mc.tag();});
                        if (count == 0)
                        {
                            mit->add_element_nonsorted(mc);
                        }
                    }
                };
                aa(m.wall_boundaries());
                aa(m.dirichlet_boundaries());
                aa(m.empty_boundaries());
                aa(m.farfield_boundaries());
            }
            else
            {
                mesh.push_front(m);
            }
        }
    }

    void SpatialPartition::set_worker_comm(const boost::mpi::communicator& comm)
    {
        worker_comm_ = comm;
    }

    void SpatialPartition::remove_dup_cells_and_points()
    {
        for (Mesh& m: mesh_)
        {
            m.remove_dup_cells_and_points();
        }
    }

    void SpatialPartition::remove_cell(const Tag& im, const Tag& ic)
    {
        mesh_p(im).remove_cell(ic);
    }

    void SpatialPartition::reset_erase_marks()
    {
        for (Mesh& m: mesh_)
        {
            m.reset_erase_marks();
        }
    }

    void SpatialPartition::mark_to_be_erased(const Tag& im, const Tag& ic)
    {
        Mesh& m = mesh_p(im);
        m.mark_to_be_erased(ic);
    }

    void SpatialPartition::erase_marked_cells()
    {
        for (Mesh& m: mesh_)
        {
            m.erase_marked_cells();
        }
    }

    void SpatialPartition::add_mesh(const SpatialPartition& other_sp)
    {
        assert(tag_ == other_sp.tag());

        for (const Mesh& om: other_sp.mesh())
        {
            auto msh = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == om.tag();});
            if (msh == mesh_.end())
            {
                add_mesh(om);
            }
        }
    }

    void SpatialPartition::add_merge_mesh_leave_dups(const SpatialPartition& other_sp)
    {
        assert(tag_ == other_sp.tag());

        for (const Mesh& om: other_sp.mesh())
        {
            auto msh = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == om.tag();});
            if (msh != mesh_.end())
            {
                //msh->merge_no_check(om);
                //std::cout << "merging - " << world_.rank() << std::endl;
                msh->merge_batch(om, world_.rank());
                //std::cout << "merged - " << world_.rank() << std::endl;
                //msh->remove_dup_cells();
                //msh->shrink_points();
                //msh->set_point_tag_index_map();
                //msh->merge_using_tags(om);
            }
            else
            {
                //std::cout << "adding - " << world_.rank() << std::endl;
                add_mesh(om);
                //std::cout << "added - " << world_.rank() << std::endl;
            }
        }
    }

    void SpatialPartition::remove_dups()
    {
        for (Mesh& msh: mesh_)
        {
            msh.remove_duplicate_cells();
            msh.remove_duplicate_points();
            msh.remove_parent_cells_of_vertices_of_all_cells();
            msh.set_parent_cell_of_vertices_of_all_cells();
            msh.update_points_from_cell_vertices();
        }
    }

    void SpatialPartition::merge_mesh(const SpatialPartition& other_sp)
    {
        assert(tag_ == other_sp.tag());

        for (const Mesh& om: other_sp.mesh())
        {
            auto msh = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == om.tag();});
            if (msh != mesh_.end())
            {
                //msh->merge_no_check(om);
                msh->merge(om);
                msh->remove_merge_duplicate_cells();
                msh->shrink_points();
                //msh->set_point_tag_index_map();
                //msh->merge_using_tags(om);
            }
        }
    }

    void SpatialPartition::merge(SpatialPartition& other_sp)
    {
        assert(tag_ == other_sp.tag());

        for (Mesh& om: other_sp.mesh_)
        {
            auto msh = std::find_if(mesh_.begin(), mesh_.end(), [&](const Mesh& m){return m.tag() == om.tag();});
            if (msh == mesh_.end())
            {
                add_mesh(om);
            }
            else
            {
                msh->remove_parent_cells_of_vertices_of_all_cells();
                msh->set_parent_cell_of_vertices_of_all_cells();
                msh->simple_merge(om);
                msh->add_points(om.cell());
                //for (const MeshCell& mc: om.cell()) {
                    //msh->add_meshpoint(mc, world_.rank());
                //}

                //msh->shrink_points();
            }
        }
    }

    /*void SpatialPartition::print_result(int iteration) const
    {
        if (master()) return;

        size_t norphan = 0.;
        for (const Mesh& m: mesh_)
        {
            for (const MeshCell& mc: m.cell())
            {
                if (mc.oga_cell_type() == OGA_cell_type_t::orphan)
                    ++norphan;
            }
        }

        std::fstream out;
        out.open("orphan-vs-iter.csv", std::fstream::out | std::fstream::app);
        if (iteration == 0)
            out << "iter,orphan\n";
        out << iteration << "," << norphan << "\n";
        out.close();

        //if (printmesh_)
        {
            int i=0;
            for (const Mesh& m: mesh_)
            {
                std::string s = "sp_";
                s.append(std::to_string(tag_()));
                s.append("_mesh_");
                s.append(std::to_string(i));
                s.append("_iter_");
                s.append(std::to_string(iteration));
                s.append(".vtk");
                m.print_as_vtk(s);
                ++i;
            }
        }
    }*/

    std::deque<std::vector<Point>> SpatialPartition::mesh_pts() const
    {
        std::deque<std::vector<Point>> pts(mesh_.size());

        int i = 0;
        for (const Mesh& m: mesh_)
        {
            for (const MeshCell& mc: m.cell())
            {
                if (!is_resident(mc.poly().centroid())) {
                    continue;
                }

                std::vector<Point>& pt = pts[i];

                for (const MeshPoint& mp: mc.point())
                {
                    pt.push_back(mp.p());
                }
            }

            ++i;
        }

        return pts;
    }

    std::deque<std::vector<Point>> SpatialPartition::mesh_pts(const Outline& outline) const
    {
        assert(false);
        std::deque<std::vector<Point>> pts(mesh_.size());

        int i = 0;
        for (const Mesh& m: mesh_)
        {
            for (const MeshCell& mc: m.cell())
            {
                //if (!is_resident(mc.poly().centroid(), outline)) {
                    //continue;
                ////////}

                std::vector<Point>& pt = pts[i];

                for (const MeshPoint& mp: mc.point())
                {
                    pt.push_back(mp.p());
                }
            }

            ++i;
        }

        return pts;
    }

    /*std::deque<AABB> SpatialPartition::mesh_aabb() const
    {
        std::deque<AABB> aabb(mesh_.size());

        for (AABB& ab: aabb)
        {
            vec3<double> minn(TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM);
            vec3<double> maxx(TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM);
            ab.set_bbox(minn, maxx);
        }

        int i = 0;
        for (const Mesh& m: mesh_)
        {
            AABB& ab = aabb[i];

            for (const MeshCell& mc: m.cell())
            {
                if (!is_resident(mc.polygon().centroid())) {
                    continue;
                }

                vec3<double> min_, max_;
                min_.set(ab.min(0), ab.min(1), ab.min(2));
                max_.set(ab.max(0), ab.max(1), ab.max(2));

                for (const MeshPoint& mp: mc.point())
                {
                    double vx = mp.p().r(0);
                    double vy = mp.p().r(1);
                    double vz = mp.p().r(2);

                    min_.set_x(std::min(min_(0), vx));
                    min_.set_y(std::min(min_(1), vy));
                    min_.set_z(std::min(min_(2), vz));
                    max_.set_x(std::max(max_(0), vx));
                    max_.set_y(std::max(max_(1), vy));
                    max_.set_z(std::max(max_(2), vz));
                }

                ab.set_bbox(min_, max_);
            }
            ++i;
        }

        return aabb;
    }*/

    /*std::deque<AABB> SpatialPartition::mesh_aabb(const Outline& outline) const
    {
        assert(false);

        std::deque<AABB> aabb(mesh_.size());

        for (AABB& ab: aabb)
        {
            vec3<double> minn(TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM, TAILOR_BIG_POS_NUM);
            vec3<double> maxx(TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM, TAILOR_BIG_NEG_NUM);
            ab.set_bbox(minn, maxx);
        }

        int i = 0;
        for (const Mesh& m: mesh_)
        {
            AABB& ab = aabb[i];

            for (const MeshCell& mc: m.cell())
            {
                //if (!is_resident(mc.polygon().centroid(), outline)) {
                    //continue;
                //}

                vec3<double> min_, max_;
                min_.set(ab.min(0), ab.min(1), ab.min(2));
                max_.set(ab.max(0), ab.max(1), ab.max(2));

                for (const MeshPoint& mp: mc.point())
                {
                    double vx = mp.p().r(0);
                    double vy = mp.p().r(1);
                    double vz = mp.p().r(2);

                    min_.set_x(std::min(min_(0), vx));
                    min_.set_y(std::min(min_(1), vy));
                    min_.set_z(std::min(min_(2), vz));
                    max_.set_x(std::max(max_(0), vx));
                    max_.set_y(std::max(max_(1), vy));
                    max_.set_z(std::max(max_(2), vz));
                }

                ab.set_bbox(min_, max_);
            }
            ++i;
        }

        return aabb;
    }*/

    void SpatialPartition::set_comm(const MPI_Comm& comm)
    {
        world_ = boost::mpi::communicator(comm, boost::mpi::comm_attach);
    }

    void SpatialPartition::set_aabb(const AABB& aabb)
    {
        aabb_ = aabb;
    }
}
