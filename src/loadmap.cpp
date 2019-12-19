#include "loadmap.h"

namespace Common
{
    void Loadmap::remove_cells_from_rm()
    {
        rm_->clear_cells();
    }

    void LM_Settings::config()
    {
        namespace po = boost::program_options;

        po::options_description all_options;

        po::options_description lm{"Load-map"};
        lm.add_options()
            ("nstripex", po::value<int>()->default_value(2), "Initial number of stripes in x-direction")
            ("nstripey", po::value<int>()->default_value(2), "Initial number of stripes in y-direction")
            ("nstripez", po::value<int>()->default_value(2), "Initial number of stripes in z-direction")
            ("refine", po::value<bool>()->default_value(true), "Can refine load-map")
            ("printlm", po::value<bool>()->default_value(false), "Print load-map")
            ("print-lm-mesh", po::value<bool>()->default_value(false), "Print meshes in load-map")
            ("refine-tol", po::value<double>()->default_value(10.), "Refinement tolerance")
            ("can-remake", po::value<bool>()->default_value(true), "Remake loadmap if degraded")
            ("force-remake", po::value<bool>()->default_value(false), "Force remake of loadmap in each time step (for developer)")
            ("refine-limit", po::value<int>()->default_value(1200), "Maximum number of refinement")
            ("balance", po::value<bool>()->default_value(true), "Refine until load balance")
            ("mergebins", po::value<bool>()->default_value(false), "Merge bins of quad-tree after refinements")
            ("load-estim", po::value<int>()->default_value(2), "Load estimation method (Default: solver)")
            ;

        all_options.add(lm);

        std::ifstream settings_file("settings.ini");

        po::store(po::parse_config_file(settings_file, all_options, true), vm_);
        po::notify(vm_);
    }

    Loadmap::Loadmap(const MPI_Comm& comm, const std::shared_ptr<Profiler>& profiler, bool verbose): profiler_(profiler), world(comm, boost::mpi::comm_attach), uniproc_(false), verbose_(verbose)
    {
        if (world.size() == 1) {
            uniproc_ = true;
        }

        read_settings();

        nworker_ = world.size() - 1;
        //reset_rm(std::ceil(std::sqrt(nworker_)));
        nrefine = 0;

        if (refine_limit_ == 0)
            dorefine = false;

        threshold_ = 10.;
    }

    void Loadmap::reset_rm(const vec3<int>& ns)
    {
        assert(ns(0) > 1);
        assert(ns(0) == ns(1));
        reset_rm(ns(0)); // there is no point of having two fields for nstripe.
    }

    void Loadmap::reset_rm(int ns)
    {
        rm_.reset();
        rm_ = std::make_unique<RegularMesh>();
        rm_->set_tag(0);
        //rm_->set_nstripe(ns, ns, ns);
        rm_->set_nstripe(ns, ns, 1);
    }

    void Loadmap::profstart(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->start(fun);
    }

    void Loadmap::profstop(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->stop(fun);
    }

    void Loadmap::recv_arrival_meshes(const SpatialPartition& sp)
    {
        if (master()) return;

        arrival_meshes.insert(arrival_meshes.end(), sp.mesh().begin(), sp.mesh().end());
    }

    void Loadmap::print(int iteration)
    {
        print_orphan(iteration);
        print_mesh(iteration);
        print_lm();
    }

    size_t Loadmap::nbin() const
    {
        return rm_->size();
    }

    void Loadmap::distribute_mti(std::vector<MeshTransferInfo>& mti) const
    {
        // prepare
        std::vector<int> mti_proc;
        std::vector<int> dest_rank(world.size());

        if (master())
        {
            //std::vector<BinRMTag> init_bin_tag;
            //init_bin_tag.resize(nworker_);

            //{
                //BinRMTag t(Tag(0), Tag(0));
                //int n = -1;
                //std::generate(init_bin_tag_.begin(), init_bin_tag_.end(),
                        //[&n, &t]()
                        //{
                        //t.set_bin_tag(Tag(++n));
                        //return t;
                        //});
            //}

            BinRMTag temp_tag(Tag(0), Tag(0));
            int pp = 0;
            //for (int p=0; p<nworker(); ++p)
            for (int p=0; p<rm_->bin().size(); ++p)
            {
                //for (const BinRMTag& t: aug_bin_tag()[p])
                {
                    //const Bin& bin = rm_->bin(t);
                    //const Bin& bin = rm_->bin(init_bin_tag[p]);
                    temp_tag.set_bin_tag(Tag(p));
                    const Bin& bin = rm_->bin(temp_tag);

                    assert(bin.parent_rm().isvalid());
                    //bin.prepare_transfer_info(p+1, mti, mti_proc, dest_rank);
                    if (p >= nworker_)
                    {
                        pp = 0;
                    }
                    bin.prepare_transfer_info(pp+1, mti, mti_proc, dest_rank, rm_->size());
                    ++pp;
                }
            }
        }

        // make sure mesh tags are unique
        for (const MeshTransferInfo& mti_: mti)
        {
            for (const MeshTransferInfo::Dest& dest: mti_.dest())
            {
                std::vector<Tag> temp;
                for (const auto& m: dest.mesh_)
                {
                    temp.push_back(m.tag_);
                }
                std::sort(temp.begin(), temp.end());
                auto it = std::adjacent_find(temp.begin(), temp.end());
                assert(it == temp.end());
            }
        }

        boost::mpi::broadcast(world, mti_proc, MASTER);
        boost::mpi::broadcast(world, dest_rank, MASTER);

        // send/recv mti to sources.
        for (int i=0; i<mti_proc.size(); ++i)
        {
            if (master())
            {
                world.send(mti_proc[i], 0, mti[i]);
            }
            else
            {
                if (world.rank() == mti_proc[i])
                {
                    MeshTransferInfo mti_temp;
                    world.recv(MASTER, 0, mti_temp);
                    mti.push_back(mti_temp);
                }
            }
        }


        std::vector<boost::mpi::request> req;
        
        // isend mti to dest.
        if (master())
        {
            {
                for (const auto& mti_: mti)
                {
                    for (int d: mti_.rank())
                    //for (const MeshTransferInfo::Dest& dest: mti_.dest())
                    {
                        //if (dest.rank_ == mti_.source()) continue;
                        if (d == mti_.source()) continue;
                        //std::cout << "I am " << mti_.source() << " sending to " << dest.rank_ << " bin tag " << dest.sp_tag_() << std::endl;
                        //req.push_back(world.isend(dest.rank_, 0, mti_));
                        req.push_back(world.isend(d, 0, mti_));
                    }
                }
            }
        }

        //world.barrier();

        // irecv mti for dest.
        if (!master())
        {
            for (int i=0; i<dest_rank.size(); ++i)
            {
                if (world.rank() == i)
                {
                    for (int j=0; j<dest_rank[i]; ++j)
                    {
                        MeshTransferInfo mti_temp;
                        //boost::mpi::status = world.recv(MASTER, 0, mti_temp);
                        world.recv(MASTER, 0, mti_temp);
                        //std::cout << "I am " << world.rank() << " receiving from " << status.source() << " bin tag " << mti_temp << std::endl;
                        //std::cout << "I am " << world.rank() << " receiving from master" << std::endl;
                        mti.push_back(mti_temp);
                    }
                }
            }
        }

        boost::mpi::wait_all(req.begin(), req.end());

        /*if (world.rank() == 1)
        {
            {
                for (const auto& mti_: mti)
                {
                    for (const MeshTransferInfo::Dest& dest: mti_.dest())
                    {
                        //if (dest.rank_ == mti_.source()) continue;
                        if (dest.rank_ != world.rank()) continue;
                        std::cout << "I am " << dest.rank_ << " receiving from " << mti_.source() << " bin tag " << dest.sp_tag_() << std::endl;
                    }
                }
            }
        }*/
        
    }

    void Loadmap::distribute_mti_2(std::vector<MeshTransferInfo>& mti) const
    {
        // prepare
        std::vector<int> mti_proc;
        std::vector<int> dest_rank(world.size());

        /*if (master())
        {
            for (int p=0; p<nworker(); ++p) {
                for (const BinRMTag& t: aug_bin_tag_[p]) {
                    std::cout << "aug_bin_tag p: " << p << " " << t.rmtag()() << " " << t.bintag()() << std::endl;
                }
            }
        }*/

        if (master())
        {
            assert(!aug_bin_tag_.empty());
            for (int p=0; p<nworker(); ++p)
            {
                //assert(!aug_bin_tag_[p].empty());
                for (const BinRMTag& t: aug_bin_tag_[p])
                {
                    const Bin& bin = rm_->bin(t);
                    //assert(!bin.cell().empty());

                    assert(bin.parent_rm().isvalid());
                    bin.prepare_transfer_info(p+1, mti, mti_proc, dest_rank, rm_->size());
                }
            }

            assert(!mti.empty());
            //assert(mti_proc.empty());
        }

        // make sure mesh tags are unique
        for (const MeshTransferInfo& mti_: mti)
        {
            for (const MeshTransferInfo::Dest& dest: mti_.dest())
            {
                std::vector<Tag> temp;
                for (const auto& m: dest.mesh_)
                {
                    temp.push_back(m.tag_);
                }
                std::sort(temp.begin(), temp.end());
                auto it = std::adjacent_find(temp.begin(), temp.end());
                assert(it == temp.end());
            }
        }

        boost::mpi::broadcast(world, mti_proc, MASTER);
        boost::mpi::broadcast(world, dest_rank, MASTER);

        // send/recv mti to sources.
        for (int i=0; i<mti_proc.size(); ++i)
        {
            if (master())
            {
                world.send(mti_proc[i], 0, mti[i]);
            }
            else
            {
                if (world.rank() == mti_proc[i])
                {
                    MeshTransferInfo mti_temp;
                    world.recv(MASTER, 0, mti_temp);
                    mti.push_back(mti_temp);
                }
            }
        }

        std::vector<boost::mpi::request> req;
        
        // isend mti to dest.
        if (master())
        {
            {
                for (const auto& mti_: mti)
                {
                    for (int i=0; i<mti_.rank().size(); ++i)
                    //for (int d: mti_.rank())
                    //for (const MeshTransferInfo::Dest& dest: mti_.dest())
                    {
                        int d = mti_.rank()[i];
                        auto dst = mti_.dest()[i];
                        //if (dest.rank_ == mti_.source()) continue;
                        if (d == mti_.source()) continue;
                        //std::cout << "I am " << mti_.source() << " sending to " << d << " ncell " << dst.ncell() << std::endl;
                        //req.push_back(world.isend(dest.rank_, 0, mti_));
                        req.push_back(world.isend(d, 0, mti_));
                    }
                }
            }
        }

        //world.barrier();

        // irecv mti for dest.
        if (!master())
        {
            for (int i=0; i<dest_rank.size(); ++i)
            {
                if (world.rank() == i)
                {
                    for (int j=0; j<dest_rank[i]; ++j)
                    {
                        MeshTransferInfo mti_temp;
                        //boost::mpi::status = world.recv(MASTER, 0, mti_temp);
                        world.recv(MASTER, 0, mti_temp);
                        //std::cout << "I am " << world.rank() << " receiving from " << status.source() << " bin tag " << mti_temp << std::endl;
                        //std::cout << "I am " << world.rank() << " receiving from master" << std::endl;
                        mti.push_back(mti_temp);
                    }
                }
            }
        }

        boost::mpi::wait_all(req.begin(), req.end());
    }

    void Loadmap::clear()
    {
        unsigned int nsx = rm_->nstripe(0);
        unsigned int nsy = rm_->nstripe(1);
        unsigned int nsz = rm_->nstripe(2);
        //donor_info_.clear();
        nrefine = 0;
        rm_.reset();
        rm_ = std::make_unique<RegularMesh>();
        //rm_ = RegularMesh();
        rm_->set_nstripe(nsx, nsy, nsz);
        rm_->set_tag(0);
        mesh_.clear();
        mesh_tag_index_map.clear();
        aug_bin_tag_.clear();
        bintag_proc_map_.clear();
    }

    void Loadmap::move_mesh(const std::vector<vec3<double>>& v)
    {
        //for (int i=0; i<mesh_.size(); ++i)
        int i = 0;
        for (Mesh& m: mesh_)
        {
            //mesh_[i].move(v[i]);
            m.move(v[i]);
            ++i;
        }
    }

    void Loadmap::rotate_mesh(const std::vector<double>& rotation, int axis, const std::vector<vec3<double>>& rot_point)
    {
        int i = 0;
        for (Mesh& m: mesh_)
        {
            //for (int i=0; i<mesh_.size(); ++i)
            //mesh_[i].rotate(rotation[i], rot_axis[i]);
            m.rotate(rotation[i], axis, rot_point[i]);
            ++i;
        }
    }

    void Loadmap::print_orphan(int iteration) const
    {
        if (!master()) return;

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
    }

    void Loadmap::print_mesh(int iteration) const
    {
        if (master()) return;

        if (printmesh_)
        {
            std::cout << "printing lm mesh" << std::endl;
            for (const Mesh& m: mesh_)
            {
                std::string s = "lm_rank_";
                s.append(std::to_string(world.rank()));
                s.append("_mesh_");
                s.append(std::to_string(m.tag()()));
                s.append("_iter_");
                s.append(std::to_string(iteration));
                s.append(".vtk");
                m.print_as_vtk(s);
            }
        }
    }

    bool Loadmap::resolution_is_enough(BinRMTag& heaviest_bt, int nmake_map)
    {
        bool ire = true;
        bool collapse = true;

        std::vector<double> sorted_load;
        std::vector<BinRMTag> sorted_tags;
        sort_load(sorted_load, sorted_tags, verbose_);

        if (master())
        {
            heaviest_bt = get_heaviest_bin(sorted_tags, verbose_);

            if (nrefine == 0)
            {
                graph_ = std::make_unique<Graph>(sorted_load[0], sorted_load[1], sorted_load[2], sorted_load[3]);
            }

            if (sorted_load.size() < nworker_) {
                collapse = false;
                ire = false;
            }

            if (balance_ == false)
            {
                if (sorted_load.size() >= nworker_)
                {
                    collapse = true;
                    ire = true;
                }
            }

            if (collapse)
            {
                graph_->connect();
                //graph_->print();
                double dev;
                bool balance = graph_->partition(nworker_, dev, aug_bin_tag_);
                assert(!aug_bin_tag_.empty());
                if (verbose_) {
                    std::cout << "dev: " << dev << std::endl;
                }

                collapse = true;

                if (balance == false) {
                    collapse = false;
                    ire = false;
                }

                std::ofstream out;
                std::string fname = "dev-vs-refine-";
                fname.append(std::to_string(nmake_map));
                fname.append(".csv");
                out.open(fname, std::fstream::app);
                if (nrefine == 0)
                    out << "ref,dev\n";

                if (collapse)
                {
                    ire = true;
                    out << nrefine << "," << dev << "\n";
                }

                out.close();

                if (dev > refine_tol_) {
                    collapse = false;
                    ire = false;
                }

                if (balance_ == false)
                {
                    std::cout << "nworker: " << nworker_ << std::endl; 
                    std::cout << "sorted_load.size(): " << sorted_load.size() << std::endl; 
                    if (sorted_load.size() >= nworker_)
                    {
                        collapse = true;
                        ire = true;
                    }
                }
            }
        }

        broadcast(world, ire, MASTER);
        broadcast(world, heaviest_bt, MASTER);

        return ire;
    }

    bool Loadmap::resolution_is_enough(BinRMTag& heaviest_bt, int nmake_map, int iteration)
    {
        bool ire;

        std::vector<double> sorted_load;
        std::vector<BinRMTag> sorted_tags;
        sort_load(sorted_load, sorted_tags, verbose_);
        if (master())
        {
            heaviest_bt = get_heaviest_bin(sorted_tags, verbose_);

            std::vector<double> aug_sorted_load;
            aug_bin_tag_.clear();
            assert(!sorted_tags.empty());
            bool collapsed = greedy_partition(nworker_, sorted_load, aug_sorted_load, sorted_tags, aug_bin_tag_, verbose_);
            double dev;
            if (!aug_sorted_load.empty()) {
                dev = estimated_deviation(aug_sorted_load);
            }
            final_dev_ = dev;

            std::ofstream out;
            std::string fname = "dev-vs-refine-";
            fname.append(std::to_string(nmake_map));
            fname.append(".csv");
            out.open(fname, std::fstream::app);
            if (nrefine == 0)
                out << "ref,dev\n";

            if (!collapsed)
            {
                ire = false;
            }
            else
            {
                ire = true;
                out << nrefine << "," << dev << "\n";
                out.close();
            }

            out.close();

            if (verbose_) {
                std::cout << "lm master outputs in reso - " << world.rank() << std::endl;
            }

            if (dev > refine_tol_)
                ire = false;
        }

        if (!uniproc_)
        {
            broadcast(world, ire, MASTER);
            broadcast(world, heaviest_bt, MASTER);
        }

        return ire;
    }

    Tag Loadmap::generate_meshtag() const
    {
        int i = BIG_NEG_NUM;
        Tag t;

        if (mesh_.empty())
        {
            t.set(0);
            return t;
        }
        for (const Mesh& m: mesh_)
        {
            i = std::max(i, m.tag()()); 
        }

        t.set(i+1);

        return t;
    }

    void Loadmap::read_settings()
    {
        LM_Settings settings;
        settings.config();

        int nsx = settings.vm_["nstripex"].as<int>();
        int nsy = settings.vm_["nstripey"].as<int>();
        int nsz = settings.vm_["nstripez"].as<int>();
        start_nstripe_.set(nsx, nsy, nsz);
        dorefine = settings.vm_["refine"].as<bool>();
        printlm_ = settings.vm_["printlm"].as<bool>();
        printmesh_ = settings.vm_["print-lm-mesh"].as<bool>();
        refine_tol_ = settings.vm_["refine-tol"].as<double>();
        refine_limit_ = settings.vm_["refine-limit"].as<int>();
        int load_estim = settings.vm_["load-estim"].as<int>();
        balance_ = settings.vm_["balance"].as<bool>();
        merge_bins_ = settings.vm_["mergebins"].as<bool>();

        if (load_estim == 0)
        {
            load_estim_type_ = LoadEstimType::area;
        }
        else if (load_estim == 1)
        {
            load_estim_type_ = LoadEstimType::hybrid;
        }
        else if (load_estim == 2)
        {
            load_estim_type_ = LoadEstimType::solver;
        }
        else if (load_estim == 3)
        {
            load_estim_type_ = LoadEstimType::minmesh;
        }
        else
        {
            assert(false);
        }
    }

    AABB Loadmap::max_aabb() const
    {
        AABB aabb;

        for (const Mesh& m: mesh_)
        {
            aabb.extend(AABB(m.rawpoint()));
        }

        //std::cout << "aabbmin0: " << aabb.min(0) << std::endl;
        //std::cout << "aabbmin1: " << aabb.min(1) << std::endl;
        //std::cout << "aabbmax0: " << aabb.max(0) << std::endl;
        //std::cout << "aabbmax1: " << aabb.max(1) << std::endl;

        return aabb;
    }

    void Loadmap::make_regular_map_uniproc_resident()
    {
        vec3<double> global_steplength_;
        vec3<double> global_aabb_min_;
        
        AABB aabb = max_aabb();

        rm_->set_aabb(aabb);
        rm_->calc_step_length();

        global_aabb_min_ = aabb.min();
        global_steplength_ = rm_->h();

        //rm_->set_aabb_min(global_aabb_min_);

        rm_->insert_bins(0, mesh_.size());

        for (const Mesh& mb: mesh_) {
            rm_->register_resident_mesh(mb, true, world.rank());
        }
    }

    void Loadmap::make_regular_map_uniproc()
    {
        vec3<double> global_steplength_;
        vec3<double> global_aabb_min_;
        
        AABB aabb = max_aabb();

        rm_->set_aabb(aabb);
        rm_->calc_step_length();

        global_aabb_min_ = aabb.min();
        global_steplength_ = rm_->h();

        //rm_->set_aabb_min(global_aabb_min_);

        rm_->insert_bins(0, mesh_.size());

        for (const Mesh& mb: mesh_) {
            rm_->register_mesh(mb, true, world.rank(), true, nullptr);
        }
    }

    void Loadmap::make_regular_map_resident()
    {
        vec3<double> global_steplength_;
        vec3<double> global_aabb_min_;
        
        std::vector<AABB> aabbs;
        AABB aabb = max_aabb();
        boost::mpi::gather(world, aabb, aabbs, MASTER);
        if (master())
        {
            for (const AABB& ab: aabbs)
            {
                aabb.extend(ab);
            }

            rm_->set_aabb(aabb);
            rm_->calc_step_length();

            global_aabb_min_ = aabb.min();
            global_steplength_ = rm_->h();

            rm_->insert_bins(0, 0);
        }

        broadcast(world, global_steplength_, MASTER);
        broadcast(world, global_aabb_min_, MASTER);

        if (!master())
        {
            rm_->set_h(global_steplength_);
            //rm_->set_aabb_min(global_aabb_min_);
            rm_->set_aabb(AABB(global_aabb_min_, vec3<double>())); // max does not matter.
            rm_->calc_aabb_max();

            rm_->insert_bins(0, mesh_.size());

            assert(!mesh_.empty());

            for (const Mesh& mb: mesh_)
            {
                assert(!mb.cell().empty());
                rm_->register_resident_mesh(mb, true, world.rank());
            }
        }
    }

    void Loadmap::make_regular_map()
    {
        if (verbose_) {
            std::cout << "making rm - " << world.rank() << std::endl;
        }

        vec3<double> global_steplength_;
        vec3<double> global_aabb_min_;
        
        std::vector<AABB> aabbs;
        AABB aabb = max_aabb();

        if (verbose_) {
            std::cout << "gathering aabss - " << world.rank() << std::endl;
        }

        boost::mpi::gather(world, aabb, aabbs, MASTER);

        if (verbose_) {
            std::cout << "gathered aabss - " << world.rank() << std::endl;
        }

        if (master())
        {
            for (const AABB& ab: aabbs)
            {
                aabb.extend(ab);
            }

            rm_->set_aabb(aabb);
            rm_->calc_step_length();

            global_aabb_min_ = aabb.min();
            global_steplength_ = rm_->h();

            if (verbose_) {
                std::cout << "inserting bins - " << world.rank() << std::endl;
            }

            rm_->insert_bins(0, 0);
        }

        broadcast(world, global_steplength_, MASTER);
        broadcast(world, global_aabb_min_, MASTER);

        if (verbose_) {
            std::cout << "broadcasted - " << world.rank() << std::endl;
        }

        if (!master())
        {
            rm_->set_h(global_steplength_);
            //rm_->set_aabb_min(global_aabb_min_);
            rm_->set_aabb(AABB(global_aabb_min_, vec3<double>())); // max does not matter.
            rm_->calc_aabb_max();

            rm_->insert_bins(0, mesh_.size());

            if (verbose_) {
                std::cout << "inserted bins - " << world.rank() << std::endl;
            }

            assert(!mesh_.empty());

            for (const Bin& b: rm_->bin())
            {
                assert(b.mesh_tag_index_map().size() == b.mesh_load().size());
            }

            for (const Mesh& mb: mesh_)
            {
                assert(!mb.cell().empty());
                rm_->register_mesh(mb, true, world.rank(), true, nullptr);
            }

            if (verbose_) {
                std::cout << "registered meshes - " << world.rank() << std::endl;
            }
        }
    }

    void Loadmap::add_mesh(const Mesh& m)
    {
        if (master())
        {
            return;
        }

        if (verbose_) {
            std::cout << "adding mesh - " << world.rank() << std::endl;
        }
        //for (const Mesh& m: mesh)
        {
            //assert(m.parent_mesh().isvalid());
            //assert(m.parent_mesh()() != -1);
            mesh_.push_back(m);
            assert(m.tag().isvalid());
            //if (m.tag()() == -1)
            //{
                //Tag t = generate_meshtag();
                //mesh_.back().set_tag(t);
            //}
            mesh_tag_index_map.insert(boost::bimap<int, int>::value_type(mesh_.back().tag()(), mesh_.size() - 1));
        }
        if (verbose_) {
            std::cout << "added mesh - " << world.rank() << std::endl;
        }
    }

    /*void Loadmap::add_mesh(const std::deque<Mesh>& mesh)
    {
        if (verbose_) {
            std::cout << "adding mesh - " << world.rank() << std::endl;
        }
        for (const Mesh& m: mesh)
        {
            if (master())
            {
                assert(false);
            }
            assert(m.parent_mesh().isvalid());
            //assert(m.parent_mesh()() != -1);
            mesh_.push_back(m);
            assert(m.tag().isvalid());
            //if (m.tag()() == -1)
            //{
                //Tag t = generate_meshtag();
                //mesh_.back().set_tag(t);
            //}
            mesh_tag_index_map.insert(boost::bimap<int, int>::value_type(mesh_.back().tag()(), mesh_.size() - 1));
        }
        if (verbose_) {
            std::cout << "added mesh - " << world.rank() << std::endl;
        }
    }*/

    /*void Loadmap::add_mesh(Mesh& m)
    {
        if (m.tag()() == -1)
        {
            Tag t = generate_meshtag();
            m.set_tag(t);
        }
        assert(m.parent_mesh()() != -1);
        mesh_.push_back(m);
        mesh_tag_index_map.insert(boost::bimap<int, int>::value_type(m.tag()(), mesh_.size() - 1));
    }*/

    void Loadmap::refine(int nmake_map, int iteration)
    {
        nrefine = 0;
        BinRMTag heaviest_bt;
        bool ire;
        if (merge_bins_)
        {
            ire = resolution_is_enough(heaviest_bt, nmake_map, iteration);
        }
        else
        {
            ire = resolution_is_enough(heaviest_bt, nmake_map);
        }

        if (dorefine)
        {
            while (!ire)
            {
                ++nrefine;
                int new_rmtag;
                int new_bt;
                rm_->refine_adaptive(mesh_, heaviest_bt, world.rank(), uniproc_, new_rmtag, new_bt);
                /*if (load_part_ == LoadPart::spatial) {
                    if (master())
                    {
                        remove_cells_from_rm();
                    }
                    gather_regular_maps();
                }*/

                if (verbose_) {
                    std::cout << "lm refined adaptive - " << world.rank() << std::endl;
                }

                if (merge_bins_)
                {
                    ire = resolution_is_enough(heaviest_bt, nmake_map, iteration);
                }
                else
                {
                    //std::vector<double> loads(4);
                    //std::vector<double> loads_r(4);
                    std::vector<double> loads;
                    std::vector<double> loads_r;
                    LoadCalculator lc(load_estim_type_, area_rep_, merge_bins_, world);
                    lc.gather_load(*rm_->bin(heaviest_bt).rm(), loads, loads_r, mesh_);

                    if (master())
                    {
                        //LoadCalculator lc(load_estim_type_, world);
                        //loads[0] = lc.load(rm_->bin(heaviest_bt).rm()->bin(0), mesh_);
                        //loads[1] = lc.load(rm_->bin(heaviest_bt).rm()->bin(1), mesh_);
                        //loads[2] = lc.load(rm_->bin(heaviest_bt).rm()->bin(2), mesh_);
                        //loads[3] = lc.load(rm_->bin(heaviest_bt).rm()->bin(3), mesh_);
                        //std::cout << "loads[0]: " << loads[0] << std::endl;
                        //std::cout << "loads[1]: " << loads[1] << std::endl;
                        //std::cout << "loads[2]: " << loads[2] << std::endl;
                        //std::cout << "loads[3]: " << loads[3] << std::endl;
                        graph_->refine(heaviest_bt.bintag()(), heaviest_bt.rmtag()(), new_rmtag, new_bt, loads[0], loads[1], loads[2], loads[3]);
                    }

                    ire = resolution_is_enough(heaviest_bt, nmake_map);
                }

                if (verbose_) {
                    std::cout << "lm reso is checked - " << world.rank() << std::endl;
                }

                if (nrefine == refine_limit_)
                    break;
            }
        }

        //if (master()) {
            //graph_->print();
        //}
        graph_.reset();
        if (master())
        {
            remove_cells_from_rm();
        }

        if (verbose_) {
            std::cout << "lm out of dorefine block - " << world.rank() << std::endl;
        }

        if (!uniproc_) {
            broadcast(world, aug_bin_tag_, MASTER);
        }

        for (const auto& a: aug_bin_tag_)
        {
            for (const auto& brmt: a)
            {
                assert(brmt.isvalid());
            }
        }

        if (master())
        {
            std::ofstream out;
            std::string fname = "refine-per-step.csv";
            out.open(fname, std::fstream::app);
            if (iteration == 0)
                out << "step,nref\n";

            out << iteration << "," << nrefine << "\n";
            out.close();
        }
    }

    void Loadmap::make_map(int iteration, int nmake_map, bool adaptive, RegType regtype)
    {
        if (verbose_) {
            std::cout << "making_lm - " << world.rank() << std::endl;
        }

        if (verbose_) {
            std::cout << "mesh size: " << mesh_.size() << " " << world.rank() << std::endl;
        }

        if (!master()) {
            assert(!mesh_.empty());
        }
        else
        {
            assert(mesh_.empty());
        }

        if (verbose_) {
            std::cout << "starting add-mesh - " << world.rank() << std::endl;
        }

        if (verbose_) {
            std::cout << "added_lm_mesh - " << world.rank() << std::endl;
        }

        if (adaptive) {
            reset_rm(start_nstripe_);
        }
        else {
            reset_rm(std::ceil(std::sqrt(nworker_)));
        }

        if (iteration == 0)
        {
            if (master())
            {
                BOOST_LOG_TRIVIAL(info) << "initial nstripe: " << rm_->nstripe(0) << "x" << rm_->nstripe(1);
                BOOST_LOG_TRIVIAL(info) << "can refine: " << dorefine;
            }
        }

        if (verbose_) {
            std::cout << "reseted rm - " << world.rank() << std::endl;
        }

        if (regtype == RegType::aabb) {
            make_regular_map();
        }
        else if (regtype == RegType::centroid) {
            make_regular_map_resident();
        }
        else {
            assert(false);
        }

        /*if (adaptive) {
            make_regular_map();
        }
        else {
            make_regular_map_resident();
        }*/
        if (verbose_) {
            std::cout << "lm made regular map - " << world.rank() << std::endl;
        }

        if (adaptive)
        {
            profstart("refine");
            refine(nmake_map, iteration);
            profstop("refine");
        }
        else
        {
            aug_bin_tag_.resize(nworker_);
            std::vector<BinRMTag> tag;
            rm_->get_bin_tags(tag);
            int j = 0;;
            for (int i=0; i<tag.size(); ++i)
            {
                if (i >= aug_bin_tag_.size()) {j = 0;}
                aug_bin_tag_[j].push_back(tag[i]);
                ++j;
            }
        }

        if (verbose_) {
            std::cout << "lm refined - " << world.rank() << std::endl;
        }
        gather_regular_maps();
        if (verbose_) {
            std::cout << "lm gather rm - " << world.rank() << std::endl;
        }
        get_bin_to_proc();
        if (verbose_) {
            std::cout << "lm getbintoproc - " << world.rank() << std::endl;
        }

        print(iteration);
    }

    /*void Loadmap::make_map_no_read(const std::deque<Mesh>& mesh, int iteration, int nmake_map, bool adaptive, RegType regtype)
    {
        if (verbose_) {
            std::cout << "making_lm - " << world.rank() << std::endl;
        }

        if (verbose_) {
            std::cout << "mesh size: " << mesh.size() << " " << world.rank() << std::endl;
        }

        if (!master()) {
            assert(!mesh.empty());
        }
        else
        {
            assert(mesh.empty());
        }

        if (verbose_) {
            std::cout << "starting add-mesh - " << world.rank() << std::endl;
        }

        profstart("add-mesh");
        add_mesh(mesh);
        profstop("add-mesh");

        if (verbose_) {
            std::cout << "added_lm_mesh - " << world.rank() << std::endl;
        }

        if (adaptive) {
            reset_rm(start_nstripe_);
        }
        else {
            reset_rm(std::ceil(std::sqrt(nworker_)));
        }

        if (iteration == 0)
        {
            if (master())
            {
                BOOST_LOG_TRIVIAL(info) << "initial nstripe: " << rm_->nstripe(0) << "x" << rm_->nstripe(1);
                BOOST_LOG_TRIVIAL(info) << "can refine: " << dorefine;
            }
        }

        if (verbose_) {
            std::cout << "reseted rm - " << world.rank() << std::endl;
        }

        if (regtype == RegType::aabb) {
            make_regular_map();
        }
        else if (regtype == RegType::centroid) {
            make_regular_map_resident();
        }
        else {
            assert(false);
        }

        //if (adaptive) {
            //make_regular_map();
        //}
        //else {
            //make_regular_map_resident();
        //}
        if (verbose_) {
            std::cout << "lm made regular map - " << world.rank() << std::endl;
        }

        if (adaptive)
        {
            profstart("refine");
            refine(nmake_map, iteration);
            profstop("refine");
        }
        else
        {
            aug_bin_tag_.resize(nworker_);
            std::vector<BinRMTag> tag;
            rm_->get_bin_tags(tag);
            int j = 0;
            for (int i=0; i<tag.size(); ++i)
            {
                if (i >= aug_bin_tag_.size()) {j = 0;}
                aug_bin_tag_[j].push_back(tag[i]);
                ++j;
            }
        }

        if (verbose_) {
            std::cout << "lm refined - " << world.rank() << std::endl;
        }
        gather_regular_maps();
        if (verbose_) {
            std::cout << "lm gather rm - " << world.rank() << std::endl;
        }
        get_bin_to_proc();
        if (verbose_) {
            std::cout << "lm getbintoproc - " << world.rank() << std::endl;
        }

        print(iteration);
    }*/

    void Loadmap::gather_regular_maps()
    {
        std::vector<RegularMesh> worker_regular_map;

        if (verbose_) {
            std::cout << "lm gathering - " << world.rank() << std::endl;
        }

        boost::mpi::gather(world, *rm_, worker_regular_map, MASTER);
        if (verbose_) {
            std::cout << "lm gathered - " << world.rank() << std::endl;
        }
        if (master())
            assert(worker_regular_map.size() == world.size());

        if (!master()) return;

        assert(!worker_regular_map.empty());

        for (const RegularMesh& _rm: worker_regular_map)
        {
            if (&_rm - worker_regular_map.data() == 0) continue;
            rm_->merge(_rm);
        }

        for (const Bin& b: rm_->bin())
        {
            if (!b.mesh_tag_index_map().left.empty())
            {
                if (b.cell().empty())
                {
                    std::cout << "bintag: " << b.tag()() << std::endl;
                    std::cout << "rmtag: " << b.parent_rm()() << std::endl;
                }
                assert(!b.cell().empty());
            }
        }
    }

    void Loadmap::sort_load(std::vector<double>& sorted_load, std::vector<BinRMTag>& sorted_tags, bool verbose)
    {
        if (verbose_) {
            std::cout << "lm calculating dev-in-load - " << world.rank() << std::endl;
        }
        verbose = false;

        LoadCalculator lc(load_estim_type_, area_rep_, merge_bins_, world);
        std::vector<double> load_global_r;
        sorted_tags = lc.sorted_bin_tags(*rm_, sorted_load, load_global_r, mesh_);

        //auto minmax = std::minmax_element(sorted_load.begin(), sorted_load.end());
        //if (*minmax.second != 0)
        //{
            //for (const auto& aa: sorted_load) {
                //std::cout << "AAAAAAAAAAAAAAAAAAAA: " << aa << std::endl;
            //}
        //}
        //assert(*minmax.second != 0);

        if (verbose_) {
            std::cout << "lm sorted tags - " << world.rank() << std::endl;
        }

        if (master())
        {
            assert(!sorted_load.empty());
            assert(!sorted_tags.empty());
        }

        if (master()) {
            std::sort(sorted_load.begin(), sorted_load.end(), std::greater<double>());
        }
    }

    BinRMTag Loadmap::get_heaviest_bin(const std::vector<BinRMTag>& sorted_tags, bool verbose)
    {
        return sorted_tags.front();

        if (verbose_) {
            std::cout << "lm sorted - " << world.rank() << std::endl;
        }
    }

    double Loadmap::estimated_deviation(const std::vector<double>& aug_sorted_load)
    {
        assert(!aug_sorted_load.empty());
        auto result = std::minmax_element(aug_sorted_load.begin(), aug_sorted_load.end());
        double minload = *result.first;
        double maxload = *result.second;
        double dev = (maxload - minload) * 100. / maxload;
        return dev;
    }

    void Loadmap::get_bin_to_proc()
    {
        if (master())
        {
            assert(!aug_bin_tag_.empty());

            for (int j=0; j<aug_bin_tag_.size(); ++j)
            {
                //std::cout << "proc - " << j+1 << std::endl;
                for(const BinRMTag& i: aug_bin_tag_[j])
                {
                    bintag_proc_map_.insert(std::pair<BinRMTag, int>(i, j+1));
                    //std::cout << i.rmtag()() << " " << i.bintag()() << std::endl;
                }
            }
        }

        if (!uniproc_) {
            broadcast(world, bintag_proc_map_, MASTER);
        }
    }

    void Loadmap::print_lm() const
    {
        if (!printlm_) {return;}
        
        std::cout << "printing lm" << std::endl;
        print_rm(*rm_);
    }

    void print_rm(const RegularMesh& rm, std::string ss, const std::map<BinRMTag, int>& bintag_proc_map_)
    {
        std::string s = ss;
        s.append(std::to_string(rm.tag()()));
        s.append(".vtk");
        rm.print(s, bintag_proc_map_);

        for (const Bin& b: rm.bin())
        {
            if (b.rm() == nullptr) continue;
            print_rm(*b.rm(), ss , bintag_proc_map_);
        }
    }

    void Loadmap::print_rm(const RegularMesh& rm) const
    {
        if (!master()) return;

        std::string s = "lm_";
        s.append(std::to_string(rm.tag()()));
        s.append(".vtk");
        rm.print(s, bintag_proc_map_);

        for (const Bin& b: rm.bin())
        {
            if (b.rm() == nullptr) continue;
            print_rm(*b.rm());
        }
    }

    bool spatial_partition(int nhead, const std::vector<double>& input, std::vector<double>& output, const std::vector<BinRMTag>& input_bin_tag, std::vector<std::vector<BinRMTag>>& output_bin_tag, bool verbose_)
    {
        if (input.size() < nhead)
            return false;

        output.reserve(nhead);
        output_bin_tag.resize(nhead);

        for (int i=0; i<nhead; ++i)
        {
            output.push_back(input[i]);
            output_bin_tag[i].push_back(input_bin_tag[i]);
        }

        if (input.size() == nhead)
            return true;

        for (auto i=input.begin()+nhead; i!=input.end(); ++i)
        {
            int dist = std::distance(input.begin(), i);
            double current = BIG_POS_NUM;
            int j = -1;
            for (int k=0; k<nhead; ++k)
            {
                if (output[k] < current)
                {
                    current = output[k];
                    j = k;
                }
            }
            if (j == -1)
            {
                if (verbose_) {
                    std::cout << "output.size = " << output.size() << std::endl;
                    std::cout << "current = " << current << std::endl;
                    std::cout << "nhead = " << nhead << std::endl;
                }
            }
            assert(j != -1);
            if (*i != 0)
            {
                output[j] += *i;
                output_bin_tag[j].push_back(input_bin_tag[dist]);
            }
        }

        int k = 0;
        for (auto i=input.rbegin(); i!=input.rend()-nhead; ++i)
        {
            if (*i != 0) continue;
            int dist = std::distance(input.begin(), i.base()) - 1;

            output[k] += *i;
            output_bin_tag[k].push_back(input_bin_tag[dist]);

            ++k;
            if (k >= nhead) k = 0;
        }

        return true;
    }

    bool greedy_partition(int nhead, const std::vector<double>& input, std::vector<double>& output, const std::vector<BinRMTag>& input_bin_tag, std::vector<std::vector<BinRMTag>>& output_bin_tag, bool verbose_)
    {
        if (input.size() < nhead)
            return false;

        output.reserve(nhead);
        output_bin_tag.resize(nhead);

        for (int i=0; i<nhead; ++i)
        {
            output.push_back(input[i]);
            output_bin_tag[i].push_back(input_bin_tag[i]);
        }

        if (input.size() == nhead)
            return true;

        for (auto i=input.begin()+nhead; i!=input.end(); ++i)
        {
            int dist = std::distance(input.begin(), i);
            double current = BIG_POS_NUM;
            int j = -1;
            for (int k=0; k<nhead; ++k)
            {
                if (output[k] < current)
                {
                    current = output[k];
                    j = k;
                }
            }
            if (j == -1)
            {
                if (verbose_) {
                    std::cout << "output.size = " << output.size() << std::endl;
                    std::cout << "current = " << current << std::endl;
                    std::cout << "nhead = " << nhead << std::endl;
                }
            }
            assert(j != -1);
            if (*i != 0)
            {
                output[j] += *i;
                output_bin_tag[j].push_back(input_bin_tag[dist]);
            }
        }

        int k = 0;
        for (auto i=input.rbegin(); i!=input.rend()-nhead; ++i)
        {
            if (*i != 0) continue;
            int dist = std::distance(input.begin(), i.base()) - 1;

            output[k] += *i;
            output_bin_tag[k].push_back(input_bin_tag[dist]);

            ++k;
            if (k >= nhead) k = 0;
        }

        return true;
    }
}
