#include "spc.h"

namespace Common
{
    void SPC_Settings::config()
    {
        namespace po = boost::program_options;

        po::options_description all_options;

        po::options_description spc{"SPC"};
        spc.add_options()
            ("print-spc-mesh", po::value<bool>()->default_value(false), "Print meshes in spc")
            ("remap-batch", po::value<bool>()->default_value(false), "Make all sp's once when sending in remapping. Batch remapping is faster but tough on memory")
            ;

        all_options.add(spc);

        std::ifstream settings_file("settings.ini");

        po::store(po::parse_config_file(settings_file, all_options, true), vm_);
        po::notify(vm_);
    }

    bool duplicate_exist(const SpatialPartitionContainer& spc)
    {
        for (const SpatialPartition& sp: spc.sp())
        {
            for (const Mesh& mesh: sp.mesh())
            {
                auto copy = mesh.cell();
                std::sort(copy.begin(), copy.end(), [&](const MeshCell& left, const MeshCell& right){return left.tag() < right.tag();});
                auto it = std::adjacent_find(copy.begin(), copy.end());
                if (it != copy.end())
                {
                    return true;
                }
            }
        }

        return false;
    }

    const std::map<BinRMTag, int>& SpatialPartitionContainer::bintag_proc_map() const
    {
        return bintag_proc_map_;
    }

    const RegularMesh& SpatialPartitionContainer::global_rm() const
    {
        return global_rm_;
    }

    void SpatialPartitionContainer::set_global_rm(const RegularMesh& rm)
    {
        global_rm_ = rm;
    }

    void SpatialPartitionContainer::make_mesh_aabbs()
    {
        for (SpatialPartition& s: sp_)
        {
            s.make_mesh_aabb();
        }
    }

    void SpatialPartitionContainer::merge_sp_meshes(std::deque<Mesh>& mesh)
    {
        if (world_.rank() != 0)
        {
            assert(!sp_.empty());
        }
        for (const SpatialPartition& s: sp_)
        {
            //s.merge_meshes_no_check(mesh);
            s.merge_meshes(mesh);
        }

        for (Mesh& m: mesh)
        {
            m.shrink_points();
            m.sort_cells();
            //m.set_cell_tag_index_map(); 
            //m.set_point_tag_index_map();
        }

        //for (Mesh& m: mesh)
        //{
            //m.remove_dup_cells_and_points();
        //}
    }

    void SpatialPartitionContainer::remap_uniproc(const std::deque<Mesh>& mesh)
    {
        SpatialPartition sp;

        AABB aabb;

        std::cout << "extending - " << world_.rank() << std::endl;

        for (const Mesh& m: mesh)
        {
            std::cout << "adding mesh to sp - " << world_.rank() << std::endl;
            sp.add_mesh(m);
            std::cout << "added mesh to sp - " << world_.rank() << std::endl;
            aabb.extend(AABB(m.rawpoint()));
            std::cout << "extented aabb - " << world_.rank() << std::endl;
        }

        std::cout << "extended - " << world_.rank() << std::endl;

        aabb.set_faces();
        aabb.set_vertices_from_bbox();
        aabb.set_faces();

        sp.set_aabb(aabb);

        add_sp(sp);
    }

    void SpatialPartitionContainer::profstart(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->start(fun);
    }

    void SpatialPartitionContainer::profstop(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->stop(fun);
    }

    void SpatialPartitionContainer::remove_dup_cells_and_points()
    {
        for (SpatialPartition& s: sp_) {
            s.remove_dup_cells_and_points();
        }
    }

    void SpatialPartitionContainer::add_sp(const std::deque<SpatialPartition>& sp, std::vector<bool>& cond)
    {
        for (int i=0; i<sp.size(); ++i)
        {
            const auto& other_sp = sp[i];
            auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});
            if (exist_sp == sp_.end())
            {
                sp_.push_back(other_sp);
                sp_.back().set_comm(world_); 
                sp_.back().set_worker_comm(worker_comm);
                sp_.back().set_profiler(profiler_);
                cond[i] = true;
                assert(!sp_.back().aabb().faces().empty());
            } 
        }
    }

    void SpatialPartitionContainer::push_sp(std::deque<SpatialPartition>& sp, const std::vector<bool>& cond)
    {
        // This function it time consuming for low number of meshes with high number of cells. Each time a cell is added to mesh, we iterate through points of mesh. So these monster meshes which usually reside in low-rank processors, cause time consumption.
        // I assumed that point insertion is the problem because we insert cells in batch with reserved space.
        // Besides point insertion, sorting of cells could also be a consumer. But sorting complexity is O(nlogn). I am not sure if partitioned sort takes same time as all in one sort. After all, n is nearly the same for all processors. 

        /*if (!master())
        {
            std::ofstream out;
            std::string fn = "incoming";
            fn.append(std::to_string(world_.rank()));
            fn.append(".csv");
            out.open(fn);
            out << sp_.size() << std::endl;
            out << sp.size() << std::endl;
            int cell_count = 0;
            int mesh_count = 0;
            for (const auto& s: sp)
            {
                mesh_count += s.mesh().size();
                for (const auto& m: s.mesh())
                {
                    cell_count += m.cell().size();
                }
            }
            out << mesh_count << std::endl;
            out << cell_count << std::endl;
            out.close();
        }*/

        for (int i=0; i<sp.size(); ++i)
        {
            if (cond[i] == true) {
                continue;
            }
            auto& other_sp = sp[i];
            auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});

            exist_sp->merge(other_sp);
            assert(!exist_sp->aabb().faces().empty());
        }
    }

    /*void SpatialPartitionContainer::make_outline(std::vector<std::vector<vec3<double>>> allpts)
    {
        Outline outline;
        if (!master())
        {
            outline.set_tag(world_.rank());
            if (!allpts.empty())
            {
                outline.build(allpts, world_.rank());
            }
        }

        boost::mpi::all_gather(world_, outline, outline_);

        if (!master())
        {
            outline_.erase(outline_.begin());
            for (Outline& ol: outline_) {
                ol.build_polygon(world_.rank());
            }
        }
    }*/

    /*void SpatialPartitionContainer::reduce_sp(std::deque<SpatialPartition>& incoming)
    {
        //std::cout << "reddd a-1 - " << world_.rank() << std::endl;
        //std::vector<vec3<double>> pts;
        std::vector<std::vector<vec3<double>>> allpts;

        auto add = [&] (const vec3<double>& v0, const vec3<double>& v1, const vec3<double>& v2, const vec3<double>& v3)
        {
            for (const auto& vvv: allpts)
            {
                if (v0 == vvv[0] && v1 == vvv[1] && v2 == vvv[2] && v3 == vvv[3])
                {
                    return;
                }
            }

            std::vector<vec3<double>> pts;
            pts.push_back(v0);
            pts.push_back(v1);
            pts.push_back(v2);
            pts.push_back(v3);
            allpts.push_back(pts);
        };

        //std::cout << "reddd a0 - " << world_.rank() << std::endl;

        if (!master())
        {
            //assert(!incoming.empty());
            if (!incoming.empty())
            {
                assert(sp_.empty());

                incoming[0].set_tag(Tag(0));
                sp_.push_back(incoming[0]);
                auto inc = incoming.begin();
                add(inc->aabb().min(), vec3<double>(inc->aabb().max(0), inc->aabb().min(1)), inc->aabb().max(), vec3<double>(inc->aabb().min(0), inc->aabb().max(1)));

                for (inc = std::next(incoming.begin()); inc != incoming.end(); ++inc)
                {
                    inc->set_tag(Tag(0));
                    assert(sp_.size() > 0);
                    sp_[0].add_merge_mesh_leave_dups(*inc);
                    sp_[0].extend_aabb(inc->aabb());
                    add(inc->aabb().min(), vec3<double>(inc->aabb().max(0), inc->aabb().min(1)), inc->aabb().max(), vec3<double>(inc->aabb().min(0), inc->aabb().max(1)));
                    //add(inc->aabb().min());
                    //add(inc->aabb().max());
                    //add(vec3<double>(inc->aabb().max(0), inc->aabb().min(1)));
                    //add(vec3<double>(inc->aabb().min(0), inc->aabb().max(1)));
                }

                sp_[0].remove_dups();
            }
        }

        //std::cout << "reddd a - " << world_.rank() << std::endl;
        //
        make_outline(allpts);
    }*/

    const std::vector<Outline>& SpatialPartitionContainer::outline() const
    {
        return outline_;
    }

    /*void SpatialPartitionContainer::reduce_sp()
    {
        if (master()) {
            return;
        }

        assert(!sp_.empty());

        SpatialPartition sp = sp_[0];

        for (int i=1; i<sp_.size(); ++i)
        {
            sp.set_tag(sp_[i].tag()); // otherwise merge_mesh() won't let the merge.
            sp.merge_mesh(sp_[i]);
            sp.extend_aabb(sp_[i].aabb());
        }

        sp_.clear();
        sp.set_tag(Tag(0));
        sp_.push_back(std::move(sp));
    }*/

    void SpatialPartitionContainer::merge_sp(const std::deque<SpatialPartition>& sp)
    {
        for (const auto& other_sp: sp)
        {
            auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});
            if (exist_sp != sp_.end())
            {
                exist_sp->merge_mesh(other_sp);
            } 
        }
    }

    /*void SpatialPartitionContainer::merge_sp(const std::deque<SpatialPartition>& sp)
    {
        for (const auto& other_sp: sp)
        {
            auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});
            if (exist_sp != sp_.end())
            {
                exist_sp->merge(other_sp);
            } 
        }
    }*/

    void SpatialPartitionContainer::add_sp(SpatialPartition& other_sp)
    {
        auto exist_sp = std::find_if(sp_.begin(), sp_.end(), [&](const auto& s){return s.tag() == other_sp.tag();});
        if (exist_sp == sp_.end())
        {
            sp_.push_back(other_sp);
            sp_.back().set_comm(world_); 
            sp_.back().set_worker_comm(worker_comm);
            sp_.back().set_profiler(profiler_);
        } 
        else
        {
            exist_sp->merge(other_sp);
        }
    }

    double SpatialPartitionContainer::dev() const
    {
        assert(master());
        return dev_;
    }

    /*void SpatialPartitionContainer::output_load(int iteration)
    {
        std::vector<double> _loads;
        gather(world_, load_, _loads, MASTER);

        if (!master()) {
            return;
        }

        double dev;
        double dev_ave;
        auto result = std::minmax_element(_loads.begin()+1, _loads.end());
        double minload = *result.first;
        double maxload = *result.second;

        double average = std::accumulate(_loads.begin()+1, _loads.end(), 0.) / (_loads.size() - 1);
        dev_ave = std::max(maxload-average, average-minload);
        dev_ave *= 100. / average;

        //dev = static_cast<double>(maxload - minload) * 100. / static_cast<double>(maxload);
        dev = (maxload - minload) * 100. / maxload;
        //for (double ll: _loads)
        //{
            //std::cout << "spc load: " << ll << std::endl;
        //}
        dev_ = dev;

        std::fstream out;
        out.open("dev-vs-iter.csv", std::fstream::out | std::fstream::app);
        if (iteration == 0)
            out << "iter,dev\n";
        out << iteration << "," << dev << "\n";
        out.close();

        out.open("devave-vs-iter.csv", std::fstream::out | std::fstream::app);
        if (iteration == 0)
            out << "iter,dev\n";
        out << iteration << "," << dev_ave << "\n";
        out.close();

        out.open("load_iter.csv", std::fstream::out | std::fstream::app);
        if (iteration == 0)
        {
            out << "iter,";
            for (int i=1; i<world_.size(); ++i)
            {
                out << "proc" << i;
                if (i != world_.size() - 1)
                    out << ",";
            }
            out << "\n";
        }
        out << iteration << ",";
        for (int i=1; i<world_.size(); ++i)
        {
            out << _loads[i];
            if (i != world_.size() - 1)
                out << ",";
        }
        out << "\n";
        out.close();
    }*/

    void SpatialPartitionContainer::read_settings()
    {
        SPC_Settings settings;
        settings.config();

        print_mesh_ = settings.vm_["print-spc-mesh"].as<bool>();
        remap_batch_ = settings.vm_["remap-batch"].as<bool>();
    }

    SpatialPartitionContainer::SpatialPartitionContainer(const MPI_Comm& comm, std::shared_ptr<Profiler> profiler, bool verbose, const std::map<BinRMTag, int>& bintag_proc_map, int nmesh): world_(comm, boost::mpi::comm_attach), profiler_(profiler), verbose_(verbose)
    {
        /*if (master())
        {
            for (const Bin& b: lm.rm().bin())
            {
                if (!b.mesh_tag_index_map().left.empty())
                {
                    assert(!b.cell().empty());
                }
            }
        }*/

        read_settings();

        bintag_proc_map_ = bintag_proc_map;
        nmesh_ = nmesh;
        boost::mpi::group world_group = world_.group();
        std::vector<int> excluded_proc = {0};
        boost::mpi::group worker_group = world_group.exclude(excluded_proc.begin(), excluded_proc.end());
        worker_comm = boost::mpi::communicator(world_, worker_group);
    }

    SpatialPartitionContainer::SpatialPartitionContainer(const MPI_Comm& comm, std::shared_ptr<Profiler> profiler, bool verbose, LoadEstimType let, int nmesh): world_(comm, boost::mpi::comm_attach), profiler_(profiler), verbose_(verbose), nmesh_(nmesh)
    {
        read_settings();

        //global_nstripe_ = lm.rm().nstripe();
        //global_steplength_ = lm.rm().h();
        //global_aabb_min_ = lm.rm().aabb().min();
        //bintag_proc_map_ = lm.bintag_proc_map();
        //nmesh_ = lm.nmesh();
        //send_proc_win = MPI_WIN_NULL;
        //send_size_win = MPI_WIN_NULL;
        //MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_2, &send_proc_win);
        //MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_wall_, &send_proc_win_wall_);
        //MPI_Win_allocate(world_.size()*sizeof(int), sizeof(int), MPI_INFO_NULL, world_, &send_proc_outer_, &send_proc_win_outer_);
        boost::mpi::group world_group = world_.group();
        std::vector<int> excluded_proc = {0};
        boost::mpi::group worker_group = world_group.exclude(excluded_proc.begin(), excluded_proc.end());
        worker_comm = boost::mpi::communicator(world_, worker_group);
    }

    void SpatialPartitionContainer::info() const
    {
        if (!master())
        {
            if (verbose_) {
                std::cout << "rank: " << world_.rank() << " | spc has " << sp_.size() << " sp." << std::endl;
            }

            for (const SpatialPartition& _sp: sp_)
            {
                _sp.info();
            }
        }
    }

    const std::vector<SpatialPartition>& SpatialPartitionContainer::sp() const
    {
        return sp_;
    }


    SpatialPartitionContainer::~SpatialPartitionContainer()
    {
        /*if (send_proc_win != MPI_WIN_NULL)
        {
            MPI_Win_free(&send_proc_win);
            MPI_Win_free(&send_proc_win_wall_);
            MPI_Win_free(&send_proc_win_outer_);
            MPI_Win_free(&send_size_win);
        }*/
    }


    void SpatialPartitionContainer::rotate_meshblocks(const Tag& _parent_mesh, double ang, int axis, const vec3<double>& rot_axis)
    {
        // this is a temporary function to mimic solver.
        // aim is to displace certain meshblocks.

        for (SpatialPartition& _sp: sp_)
        {
            _sp.rotate_mesh(_parent_mesh, ang, axis, rot_axis);
        }
    }

    void SpatialPartitionContainer::move_meshblocks(const Tag& _parent_mesh, const vec3<double>& v)
    {
        // this is a temporary function to mimic solver.
        // aim is to displace certain meshblocks.

        for (SpatialPartition& _sp: sp_)
        {
            _sp.move_mesh(_parent_mesh, v);
        }
    }

    void SpatialPartitionContainer::connect_cells_uniproc()
    {
        for (SpatialPartition& _sp: sp_) {
            _sp.connect_cells_uniproc();
        }
    }

    void SpatialPartitionContainer::connect_cells()
    {
        if (master()) return;

        for (SpatialPartition& _sp: sp_) {
            _sp.connect_cells();
        }
    }

    bool SpatialPartitionContainer::master() const 
    {
        if (world_.rank() == MASTER)
        {
            return true;
        }

        return false;
    }

    void SpatialPartitionContainer::make_regular_maps_uniproc()
    {
        for (SpatialPartition& _sp: sp_)
        {
            _sp.make_regular_maps_for_mesh_uniproc(profiler_);
        }
    }

    void SpatialPartitionContainer::make_regular_maps()
    {
        for (SpatialPartition& _sp: sp_)
        {
            _sp.make_regular_maps_for_mesh();
            //for (const auto& _rm: _sp.rm())
            //{
                //std::cout << "rank: " << world_.rank() << " rm size: " << _rm.bin().size() << std::endl;
            //}
        }
    }

    void SpatialPartitionContainer::print_all_meshes_in_partitions_uniproc()
    {
        int nmesh = 0;
        for (SpatialPartition& _sp: sp_)
        {
            nmesh += _sp.mesh().size();
            for (int i=0; i<_sp.mesh().size(); ++i)
            {
                std::string s = "spc_rank_";
                s.append(std::to_string(world_.rank()));
                s.append("_sp_");
                s.append(std::to_string(&_sp - sp_.data()));
                s.append("_mb_");
                //s.append(std::to_string(i));
                s.append(std::to_string(_sp.mesh()[i].tag()()));
                s.append(".vtk");
                _sp.mesh()[i].print_as_vtk(s);
            }
        }
    }

    void SpatialPartitionContainer::print_all_meshes_in_partitions()
    {
        if (master()) return;

        int nmesh = 0;
        for (SpatialPartition& _sp: sp_)
        {
            nmesh += _sp.mesh().size();
            for (int i=0; i<_sp.mesh().size(); ++i)
            {
                std::string s = "spc_rank_";
                s.append(std::to_string(world_.rank()));
                s.append("_sp_");
                s.append(std::to_string(&_sp - sp_.data()));
                s.append("_mb_");
                //s.append(std::to_string(i));
                s.append(std::to_string(_sp.mesh()[i].tag()()));
                s.append(".vtk");
                _sp.mesh()[i].print_as_vtk(s);
            }
        }
    }
}
