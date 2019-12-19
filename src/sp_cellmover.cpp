#include "sp_cell_mover.h"

namespace Common
{
    void SpCellMover::add_arrival_info(size_t size, int disp_begin, int disp_end)
    {
        arrival_info_.add_size(size);
        arrival_info_.add_disp(disp_begin, disp_end);
    }

    const Tag& SpCellMover::tag() const
    {
        return tag_;
    }

    void SpCellMover::push_moving_elements()
    {
        if (master()) return;

        //for (const Mesh& m: mesh_)
        //{
            //profstart("movcel-real");
            movcel_.add_cell_real(sp_->mesh_, "interior");
            wall_cell_.add_cell_real(sp_->mesh_, "wall");
            dirichlet_cell_.add_cell_real(sp_->mesh_, "dirichlet");
            empty_cell_.add_cell_real(sp_->mesh_, "empty");
            farfield_cell_.add_cell_real(sp_->mesh_, "farfield");
            //profstop("movcel-real");
        //}
    }

    bool SpCellMover::master() const
    {
        if (world_.rank() == 0) {
            return true;
        }

        return false;
    }

    void SpCellMover::clear_arrival_cell()
    {
        arrival_cell_.clear();
    }

    void SpCellMover::clear_arrival_wall_cell()
    {
        arrival_wall_cell_.clear();
    }

    void SpCellMover::clear_arrival_empty_cell()
    {
        arrival_empty_cell_.clear();
    }

    void SpCellMover::clear_arrival_farfield_cell()
    {
        arrival_farfield_cell_.clear();
    }

    void SpCellMover::clear_arrival_dirichlet_cell()
    {
        arrival_dirichlet_cell_.clear();
    }

    std::vector<MeshCell>& SpCellMover::arrival_cell()
    {
        return arrival_cell_;
    }

    std::vector<MeshCell>& SpCellMover::arrival_wall_cell()
    {
        return arrival_wall_cell_;
    }

    std::vector<MeshCell>& SpCellMover::arrival_dirichlet_cell()
    {
        return arrival_dirichlet_cell_;
    }

    std::vector<MeshCell>& SpCellMover::arrival_empty_cell()
    {
        return arrival_empty_cell_;
    }

    std::vector<MeshCell>& SpCellMover::arrival_farfield_cell()
    {
        return arrival_farfield_cell_;
    }

    SpCellMover::SpCellMover(SpatialPartition* sp, const MPI_Comm& comm, bool mergebins): world_(comm, boost::mpi::comm_attach), mergebins_(mergebins)
    {
        sp_ = sp;
        assert(sp_ != nullptr);
        tag_ = sp_->tag();
        assert(tag_.isvalid());
    }
    const MovingCell& SpCellMover::movcel() const
    {
        return movcel_;
    }

    /*void SpCellMover::reserve_moving_elements(const std::vector<Outline>& outlines)
    {
        arrival_cell_.clear();
        movcel_.clear();
        wall_cell_.clear();
        outer_cell_.clear();

        movcel_.reserve_for_mcp(sp_->mesh_.size());
        assert(!outlines.empty());

        for (const Mesh& m: sp_->mesh_)
        {
            //assert(sp_->rm(m.tag()).aabb().min(0) < sp_->rm(m.tag()).aabb().max(0));
            AABB mesh_aabb(m.rawpoint());
            auto mesh_aabb_vertices = mesh_aabb.vertices();

            auto myoutline = std::find_if(outlines.begin(), outlines.end(), [&](const Outline& ol){return ol.tag()() == world_.rank();});
            assert(myoutline != outlines.end());
            auto myoutline_aabb = AABB(*myoutline);

            // if cell is inside its current outline. 
            bool all_in = true;
            for (const auto& v: mesh_aabb_vertices)
            {
                if (!myoutline->do_contain(v, true))
                {
                    all_in = false;
                    break;
                }
            }

            if (all_in) {
                continue;
            }

            std::vector<const Outline*> potential_outline;
            for (const Outline& outline: outlines)
            {
                if (outline.tag() == myoutline->tag()) {
                    continue;
                }

                if (AABB(outline).do_intersect(mesh_aabb) == false) {
                    continue;
                }

                if (CGAL::do_intersect(mesh_aabb.cgal_polygon(), outline.cgal_polygon()))
                {
                    potential_outline.push_back(&outline);
                }
            }

            if (potential_outline.empty()) {
                continue;
            }

            //auto pol0 = m.cell().front().polygon().cgalpolygon();
            //auto pol1 = m.cell().back().polygon().cgalpolygon();

            std::vector<bool> cell_in(m.cell().size(), false);
            int cic = -1;
            for (const MeshCell& c: m.cell())
            {
                ++cic;
                AABB cell_aabb(c.polygon());
                if (!cell_aabb.do_intersect(myoutline_aabb)) {
                    continue;
                }

                auto cell_aabb_vertices = cell_aabb.vertices();
                bool cell_all_in = true;
                for (const auto& v: cell_aabb_vertices)
                {
                    if (!cell_aabb.do_intersect(myoutline_aabb))
                    {
                        cell_all_in = false;
                        break;
                    }
                }

                if (cell_all_in) {
                    cell_in[cic] = true;
                    continue;
                }

                cell_all_in = true;
                for (const auto& v: cell_aabb_vertices)
                {
                    if (!myoutline->do_contain(v, true))
                    {
                        cell_all_in = false;
                        break;
                    }
                }

                if (cell_all_in) {
                    cell_in[cic] = true;
                    continue;
                }
            }

            //for (const Outline* outline: potential_outline)
            {
                //auto outline_aabb = AABB(*outline);
                //for (const Bin& bin: sp_->rm(m.tag()).bin())
                {
                    //auto bin_vertices = bin.aabb().vertices();
                    //for (const auto& v: bin_vertices)
                    {
                        //if (outline->do_contain(v, false))
                        {
                            cic = -1;
                            for (const MeshCell& c: m.cell())
                                //for (const auto& bincell: bin.cell())
                            {
                                ++cic;
                                if (cell_in[cic]) {
                                    continue;
                                }
                                AABB cell_aabb(c.polygon());
                                auto vertices = cell_aabb.vertices();
                                for (const Outline* outline: potential_outline)
                                {
                                    auto outline_aabb = AABB(*outline);
                                    if (!cell_aabb.do_intersect(outline_aabb)) {
                                        continue;
                                    }
                                    // make dummy cgal polygons to check if construction of cgal polygons from cell_aabb and outline is costly.
                                    //CGAL::do_intersect(pol0, pol1);
                                    if (!CGAL::do_intersect(cell_aabb.cgal_polygon(), outline->cgal_polygon())) {
                                        continue;
                                    }
                                    //for (const auto& v: vertices)
                                    {
                                        //if (outline->do_contain(v, false))
                                        {
                                            //assert(m.tag()() != 0);

                                            int dest_rank = outline->tag()();
                                            assert(dest_rank != world_.rank());
                                            movcel_.add_cell_pre(dest_rank, c, dest_rank);

                                            for (const Tag& bc: c.wall_boundary())
                                            {
                                                assert(m.wall_boundary(bc).is_wall());
                                                wall_cell_.add_cell_pre(dest_rank, m.wall_boundary(bc), dest_rank);
                                            }
                                            for (const Tag& bc: c.outer_boundary())
                                            {
                                                assert(m.outer_boundary(bc).is_outer());
                                                outer_cell_.add_cell_pre(dest_rank, m.outer_boundary(bc), dest_rank);
                                            }

                                            //break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }*/

    void SpCellMover::reserve_moving_elements(const std::map<BinRMTag, int>& bintag_proc_map, const RegularMesh& global_rm_, size_t& ut_size, size_t& move_size, int& nlevel, size_t& preresi, size_t& postresi)
    {
        // This function is time consuming for sp's which have higher number of meshes. It seems like in add_cell_pre first we find_if through mcp (meshes) and if no such mcp exists we push_back. Both find_if and push_back might be the reason of time consumption.
        // I can reserve mesh_.size() capacity for mcp to avoid reallocation. I tried but no visible change in time consumption.
        // I don't think we will benefit from replacing linear search with sorted binary search because number of meshes for test_5 was maximum 26 and minimum 1.
        // Results indicate that number of meshes or sp is not the reason.
        // The reason might be get_bintag_adaptive(). Some cells fall into highly refined region of adaptive map.

        move_size = 0;
        nlevel = 0;
        preresi = 0;
        postresi = 0;

        arrival_cell_.clear();
        movcel_.clear();
        wall_cell_.clear();
        dirichlet_cell_.clear();
        empty_cell_.clear();
        farfield_cell_.clear();

        movcel_.reserve_for_mcp(sp_->mesh_.size());

        // determine displaced cells.
        for (const Mesh& m: sp_->mesh_)
        {
            for (const MeshCell& c: m.cell())
            {
                // ghosts will be added implicitly.
                //if (c.is_ghost())
                //continue;

                // we should send only residents!!!

                for (const MeshPoint& p: c.point())
                    assert(p.parent_cell().size() > 0);

                std::vector<BinRMTag> tag;
                assert(global_rm_.tag()() == 0);
                //profstart("gbta");

                ++preresi;

                if (sp_->fully_resident(c)) {
                    continue;
                }

                ++postresi;

                int temp;

                global_rm_.get_bintag_adaptive(AABB(c.poly()), tag, temp);
                nlevel += temp;
                //profstop("gbta");

                //profstart("gbta-2");
                // get rid of duplicate indices.
                std::vector<BinRMTag> ut = tag;
                std::sort(ut.begin(), ut.end(), &BinRMTag::sorter);
                ut.erase(std::unique(ut.begin(), ut.end()), ut.end());

                // determine corner bins.
                std::vector<RegularMeshIndex> minb(ut.size());
                std::vector<RegularMeshIndex> maxb(ut.size());
                for (int i=0; i<ut.size(); ++i)
                {
                    int mini = BIG_POS_NUM;
                    int minj = BIG_POS_NUM;
                    int mink = BIG_POS_NUM;
                    int maxi = BIG_NEG_NUM;
                    int maxj = BIG_NEG_NUM;
                    int maxk = BIG_NEG_NUM;
                    for (const BinRMTag& _tag: tag)
                    {
                        if (_tag.rmtag() == ut[i].rmtag())
                        {
                            const Bin& _bin = global_rm_.bin(_tag);

                            mini = std::min(mini, _bin.index().i());
                            minj = std::min(minj, _bin.index().j());
                            mink = std::min(mink, _bin.index().k());
                            maxi = std::max(maxi, _bin.index().i());
                            maxj = std::max(maxj, _bin.index().j());
                            maxk = std::max(maxk, _bin.index().k());
                        }
                    }

                    minb[i].set(mini, minj, mink);
                    maxb[i].set(maxi, maxj, maxk);
                }
                //profstop("gbta-2");

                //profstart("gbta-3");
                // loop through indices.
                ut_size = ut.size();
                for (int u=0; u<ut.size(); ++u)
                {
                    auto it = bintag_proc_map.find(ut[u]);
                    assert(it != bintag_proc_map.end());
                    unsigned int dest_rank = it->second;
                    if (dest_rank == world_.rank()) {
                        continue;
                    }

                    //for (int k=minb[u].k(); k<=maxb[u].k(); ++k)
                    {
                        //for (int i=minb[u].i(); i<=maxb[u].i(); ++i)
                        {
                            //for (int j=minb[u].j(); j<=maxb[u].j(); ++j)
                            {
                                const Bin& _bin = global_rm_.bin(ut[u]);
                                if (_bin.aabb().do_intersect(AABB(c.poly())))
                                {
                                    if (sp_->tag() != _bin.tag())
                                    {
                                        ++move_size;
                                        //auto it = bintag_proc_map.find(ut[u]);
                                        //assert(it != bintag_proc_map.end());
                                        //unsigned int dest_rank = it->second;
                                        //movcel_.add_cell(dest_rank, c, _bin.tag()());
                                        if (!mergebins_) {
                                            movcel_.add_cell_pre(dest_rank, c, _bin.tag()());
                                        }
                                        else {
                                            assert(false);
                                            //movcel_.add_cell_pre(dest_rank, c, 0); // since now we have only one sp and tags (0) must match.
                                        }
                                        //movcel_.add_cell_pre(dest_rank, c, pairing_function(tag_(), ut[u].bintag()()));

                                        // add boundaries if any.
                                        for (const Tag& bc: c.wall_boundary())
                                        {
                                            assert(m.wall_boundary(bc).boutype() == "wall");
                                            if (!mergebins_) {
                                                wall_cell_.add_cell_pre(dest_rank, m.wall_boundary(bc), _bin.tag()());
                                            }
                                            else {
                                                assert(false);
                                            }
                                        }
                                        for (const Tag& bc: c.dirichlet_boundary())
                                        {
                                            assert(m.dirichlet_boundary(bc).boutype() == "dirichlet");
                                            if (!mergebins_) {
                                                dirichlet_cell_.add_cell_pre(dest_rank, m.dirichlet_boundary(bc), _bin.tag()());
                                            }
                                            else {
                                                assert(false);
                                            }
                                        }
                                        for (const Tag& bc: c.empty_boundary())
                                        {
                                            assert(m.empty_boundary(bc).boutype() == "empty");
                                            if (!mergebins_) {
                                                empty_cell_.add_cell_pre(dest_rank, m.empty_boundary(bc), _bin.tag()());
                                            }
                                            else {
                                                assert(false);
                                            }
                                        }
                                        for (const Tag& bc: c.farfield_boundary())
                                        {
                                            assert(m.farfield_boundary(bc).boutype() == "farfield");
                                            if (!mergebins_) {
                                                farfield_cell_.add_cell_pre(dest_rank, m.farfield_boundary(bc), _bin.tag()());
                                            }
                                            else {
                                                assert(false);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                }
                }
                //profstop("gbta-3");
            }
        }
    }

    void SpCellMover::inform_receivers(MPI_Win& send_proc_win, int* send_proc_2, MPI_Win& send_size_win, int* send_size, std::vector<int>& receive_count, std::vector<int> receive_size)
    {
        if (!master())
        {
            for (int i=0; i<movcel_.cell().size(); ++i)
            {
                for (const auto a: movcel_.cell(i))
                {
                    assert(a.boutype() == "interior");
                }
                assert(!movcel_.cell(i).empty());

                /*double limit = 1e9;
                int c;
                if (movcel_.index_to_rank(i) != world_.rank())
                {
                    int count = 0;
                    double size = movcel_.cell(i).size();
                    c = std::ceil(size / limit);
                    int current = int(size);
                    int oldm = 0;

                    for (int j=0; j<c; ++j)
                    {
                        int m = std::ceil(current / (c-j));
                        current -= m;

                        oldm = m;
                    }
                }*/

                assert(movcel_.index_to_rank(i) != world_.rank());
                //if (movcel_.index_to_rank(i) != world_.rank())
                {
                    ++receive_count[movcel_.index_to_rank(i)];
                    receive_size[movcel_.index_to_rank(i)] += movcel_.cell(i).size();

                    //int a = c;
                    //int a = 1;
                    //std::cout << "insidebef informing rece - " << world_.rank() << " " << movcel_.index_to_rank(i) << " " << a << std::endl;
                    //assert(movcel_.index_to_rank(i) > 0);
                    //assert(movcel_.index_to_rank(i) < world_.size());
                    //MPI_Win_lock(MPI_LOCK_EXCLUSIVE, movcel_.index_to_rank(i), 0, send_proc_win);
                    //std::cout << "lock informing rece - " << world_.rank() << std::endl;
                    //MPI_Accumulate(&a, 1, MPI_INT, movcel_.index_to_rank(i), movcel_.index_to_rank(i), 1, MPI_INT, MPI_SUM, send_proc_win);
                    //std::cout << "unlock informing rece - " << world_.rank() << std::endl;
                    //MPI_Win_unlock(movcel_.index_to_rank(i), send_proc_win);
                    //std::cout << "insideafter informing rece - " << world_.rank() << std::endl;
//
                    //a = movcel_.cell(i).size();
                    //MPI_Win_lock(MPI_LOCK_EXCLUSIVE, movcel_.index_to_rank(i), 0, send_size_win);
                    //MPI_Accumulate(&a, 1, MPI_INT, movcel_.index_to_rank(i), movcel_.index_to_rank(i), 1, MPI_INT, MPI_SUM, send_size_win);
                    //MPI_Win_unlock(movcel_.index_to_rank(i), send_size_win);
                    //std::cout << "insideafterrr informing rece - " << world_.rank() << std::endl;
                }
                //else
                //{
                    //send_size[movcel_.index_to_rank(i)] += movcel_.cell(i).size();
                //}

                if (!mergebins_) {
                    assert(tag_() != movcel_.index_to_binindex(i)); // since now only one sp.
                }
            }
        }
    }

    void SpCellMover::try_to_receive(int* send_proc_2, std::vector<MeshCell>& arrivals_)
    {
        for (int j=0; j<send_proc_2[world_.rank()]; ++j)
        {
            if (rece[j] == false)
            {
                boost::optional<boost::mpi::status> status = world_.iprobe();
                if (status)
                {
                    if (tag_() == status->tag()) //uncomment for greedy
                    {
                        std::vector<MeshCell> temp;
                        //std::cout << world_.rank() << " is receiving from " << status->source() << " tag: " << status->tag() << " send_proc: " << send_proc_2[world_.size()] << std::endl;

                        //world_.recv(status->source(), tag_(), temp);
                        world_.recv(status->source(), status->tag(), temp);

                        //std::cout << world_.rank() << " received from " << status->source() << " with size " << temp.size() << " current arruvak size: " << arrivals_.size() << std::endl;
                        assert(!temp.empty());

                        int tempsize = temp.size();
                        arrivals_.insert(arrivals_.end(), std::make_move_iterator(temp.begin()), std::make_move_iterator(temp.end()));

                        //std::cout << world_.rank() << "'s arrival size: " << arrivals_.size() << std::endl;

                        assert(tempsize != 0);
                        add_arrival_info(tempsize, arrivals_.size() - tempsize, arrivals_.size());
                        rece[j] = true;
                        break;
                    }
                }
            }
        }
    }

    void SpCellMover::complete_receive(int* send_proc_2, std::vector<MeshCell>& arrivals_)
    {
        bool not_complete = true;
        while (not_complete)
        {
            for (int j=0; j<send_proc_2[world_.rank()]; ++j)
            {
                if (rece[j] == false)
                {
                    //std::cout << world_.rank() << " is probing for and send_proc_suze is " << send_proc_2[world_.rank()] << std::endl;
                    boost::optional<boost::mpi::status> status = world_.probe();
                    if (status)
                    {
                        if (tag_() == status->tag()) //uncomment for greedy. does greedy use this function?
                        {
                            std::vector<MeshCell> temp;
                            //std::cout << world_.rank() << " is receiving from " << status.source() << std::endl;

                            world_.recv(status->source(), status->tag(), temp);

                            //std::cout << world_.rank() << " receivedd from " << status.source() << " with size " << temp.size() << " current arruvak size: " << arrivals_.size() << std::endl;

                            assert(!temp.empty());
                            int tempsize = temp.size();
                            arrivals_.insert(arrivals_.end(), std::make_move_iterator(temp.begin()), std::make_move_iterator(temp.end()));

                            //std::cout << world_.rank() << "'s arrival size: " << arrivals_.size() << std::endl;

                            assert(tempsize != 0);
                            add_arrival_info(tempsize, arrivals_.size() - tempsize, arrivals_.size());
                            rece[j] = true;
                            break;
                        }
                    }
                }
            }

            not_complete = false;
            if (std::any_of(rece.begin(), rece.end(), [](bool b){return b == false;}))
            {
                not_complete = true;
            }
        }
    }

    void SpCellMover::send_moving_cells(int* send_proc_2, std::vector<MeshCell>& arrivals_)
    {
        rece.resize(send_proc_2[world_.rank()], false);

        for (int i=0; i<movcel_.cell().size(); ++i)
        {
            for (const auto a: movcel_.cell(i))
            {
                assert(a.boutype() == "interior");
            }
            assert(!movcel_.cell(i).empty());

            //double limit = 1e9;
            //int c;
            if (movcel_.index_to_rank(i) != world_.rank())
            {
                /*int count = 0;
                double size = movcel_.cell(i).size();
                c = std::ceil(size / limit);
                int current = int(size);
                int oldm = 0;*/

                //for (int j=0; j<c; ++j)
                {
                    //int m = std::ceil(current / (c-j));
                    //current -= m;

                    //std::cout << world_.rank() << " is sending to " << movcel_.index_to_rank(i) << " " << movcel_.index_to_binindex(i) << " " << movcel_.cell(i).size() << std::endl;
                    //std::cout << world_.rank() << " is sending to " << movcel_.index_to_rank(i) << " " << world_.rank() << " " << size << " " << c << " " << m << " " << oldm << std::endl;
                    //boost::mpi::request req = world_.isend(movcel_.index_to_rank(i), movcel_.index_to_binindex(i), &(movcel_.cell(i)[oldm]), m);
                    //boost::mpi::request req = world_.isend(movcel_.index_to_rank(i), j, &(movcel_.cell(i)[oldm]), m);
                    boost::mpi::request req = world_.isend(movcel_.index_to_rank(i), movcel_.index_to_binindex(i), movcel_.cell(i));
                    //boost::mpi::request req = world_.isend(movcel_.index_to_rank(i), world_.rank(), movcel_.cell(i));
                    //boost::mpi::request req = world_.isend(movcel_.index_to_rank(i), world_.rank(), movcel_.cell(i));
                    //std::cout << world_.rank() << " isend to " << movcel_.index_to_rank(i) << " " << movcel_.index_to_binindex(i) << " " << size << " " << c << " " << m << " " << oldm << std::endl;
                    //std::cout << world_.rank() << " isend to " << movcel_.index_to_rank(i) << " " << world_.rank() << " " << size << " " << c << " " << m << " " << oldm << std::endl;
                    //std::cout << world_.rank() << " isend to " << movcel_.index_to_rank(i) << " " << movcel_.index_to_binindex(i) << " " << movcel_.cell(i).size() << std::endl;
                    //oldm = m;
                    //std::cout << world_.rank() << " is trying to receive" << std::endl;
                    //std::cout << world_.rank() << " starting testing " << i << " " << movcel_.cell().size() << std::endl;
                    while (true)
                    {
                        if (req.test()) {
                            break;
                        }
                        try_to_receive(send_proc_2, arrivals_);
                    }
                    //std::cout << world_.rank() << " tried to receive " << i << " " << movcel_.cell().size() << std::endl;
                }
            }
        }

        //std::cout << "rank: " << world_.rank() << " completing receival" << std::endl;

        //std::cout << world_.rank() << " is completing receive" << std::endl;
        complete_receive(send_proc_2, arrivals_);
        //std::cout << world_.rank() << " completed receive" << std::endl;

        //std::cout << "rank: " << world_.rank() << " completed receival" << std::endl;
    }

    void SpCellMover::send_moving_cells(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win, int* send_proc_2, MPI_Win& send_size_win, int* send_size)
    {
        for (int i=0; i<movcel_.cell().size(); ++i)
        {
            for (const auto a: movcel_.cell(i))
            {
                assert(a.boutype() == "interior");
            }
            assert(!movcel_.cell(i).empty());

            double limit = 1e9;
            int c;
            if (movcel_.index_to_rank(i) != world_.rank())
            {
                int count = 0;
                double size = movcel_.cell(i).size();
                c = std::ceil(size / limit);
                int current = int(size);
                int oldm = 0;

                for (int j=0; j<c; ++j)
                {
                    int m = std::ceil(current / (c-j));
                    current -= m;

                    std::cout << world_.rank() << " is sending to " << movcel_.index_to_rank(i) << " " << movcel_.index_to_binindex(i) << " " << size << " " << c << " " << m << " " << oldm << std::endl;
                    request.push_back(world_.isend(movcel_.index_to_rank(i), movcel_.index_to_binindex(i), &(movcel_.cell(i)[oldm]), m));
                    oldm = m;
                }

                //request.push_back(world_.isend(movcel_.index_to_rank(i), movcel_.index_to_binindex(i), movcel_.cell(i)));
                //std::cout << world_.rank() << "sent to " << movcel_.index_to_rank(i) << " " << movcel_.index_to_binindex(i) << " " << movcel_.cell(i).size() << std::endl;
            }

            if (movcel_.index_to_rank(i) != world_.rank())
            {
                //int a = 1;
                int a = c;
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, movcel_.index_to_rank(i), 0, send_proc_win);
                MPI_Accumulate(&a, 1, MPI_INT, movcel_.index_to_rank(i), movcel_.index_to_rank(i), 1, MPI_INT, MPI_SUM, send_proc_win);
                MPI_Win_unlock(movcel_.index_to_rank(i), send_proc_win);

                a = movcel_.cell(i).size();
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, movcel_.index_to_rank(i), 0, send_size_win);
                MPI_Accumulate(&a, 1, MPI_INT, movcel_.index_to_rank(i), movcel_.index_to_rank(i), 1, MPI_INT, MPI_SUM, send_size_win);
                MPI_Win_unlock(movcel_.index_to_rank(i), send_size_win);
            }
            else
            {
                //++send_proc_2[movcel_.index_to_rank(i)];
                send_size[movcel_.index_to_rank(i)] += movcel_.cell(i).size();
            }

            if (!mergebins_) {
                assert(tag_() != movcel_.index_to_binindex(i)); // since now only one sp.
            }
        }
    }

    void SpCellMover::send_moving_farfield(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win_farfield_, int* send_proc_farfield_, std::vector<int>& receive_count)
    {
        for (int i=0; i<farfield_cell_.cell().size(); ++i)
        {
            for (const auto a: farfield_cell_.cell(i))
            {
                assert(a.boutype() == "farfield");
            }
            if (farfield_cell_.index_to_rank(i) != world_.rank()) {
                request.push_back(world_.isend(farfield_cell_.index_to_rank(i), farfield_cell_.index_to_binindex(i), farfield_cell_.cell(i)));
            }

            //++send_proc[farfield_cell_.index_to_rank(i)];
            assert(farfield_cell_.index_to_rank(i) != world_.rank());
            //if (farfield_cell_.index_to_rank(i) != world_.rank())
            {
                ++receive_count[farfield_cell_.index_to_rank(i)];
                //int a = 1;
                //MPI_Win_lock(MPI_LOCK_EXCLUSIVE, farfield_cell_.index_to_rank(i), 0, send_proc_win_farfield_);
                //MPI_Accumulate(&a, 1, MPI_INT, farfield_cell_.index_to_rank(i), farfield_cell_.index_to_rank(i), 1, MPI_INT, MPI_SUM, send_proc_win_farfield_);
                //MPI_Win_unlock(farfield_cell_.index_to_rank(i), send_proc_win_farfield_);
            }
            //else
                //++send_proc_farfield_[farfield_cell_.index_to_rank(i)];

            assert(tag_() != farfield_cell_.index_to_binindex(i));
        }
    }

    void SpCellMover::send_moving_empty(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win_empty_, int* send_proc_empty_, std::vector<int>& receive_count)
    {
        for (int i=0; i<empty_cell_.cell().size(); ++i)
        {
            for (const auto a: empty_cell_.cell(i))
            {
                assert(a.boutype() == "empty");
            }
            if (empty_cell_.index_to_rank(i) != world_.rank()) {
                request.push_back(world_.isend(empty_cell_.index_to_rank(i), empty_cell_.index_to_binindex(i), empty_cell_.cell(i)));
            }

            //++send_proc[empty_cell_.index_to_rank(i)];
            assert(empty_cell_.index_to_rank(i) != world_.rank());
            //if (empty_cell_.index_to_rank(i) != world_.rank())
            {
                ++receive_count[empty_cell_.index_to_rank(i)];
                //int a = 1;
                //MPI_Win_lock(MPI_LOCK_EXCLUSIVE, empty_cell_.index_to_rank(i), 0, send_proc_win_empty_);
                //MPI_Accumulate(&a, 1, MPI_INT, empty_cell_.index_to_rank(i), empty_cell_.index_to_rank(i), 1, MPI_INT, MPI_SUM, send_proc_win_empty_);
                //MPI_Win_unlock(empty_cell_.index_to_rank(i), send_proc_win_empty_);
            }
            //else
                //++send_proc_empty_[empty_cell_.index_to_rank(i)];

            assert(tag_() != empty_cell_.index_to_binindex(i));
        }
    }

    void SpCellMover::send_moving_dirichlet(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win_dirichlet_, int* send_proc_dirichlet_, std::vector<int>& receive_count)
    {
        for (int i=0; i<dirichlet_cell_.cell().size(); ++i)
        {
            for (const auto a: dirichlet_cell_.cell(i))
            {
                assert(a.boutype() == "dirichlet");
            }
            if (dirichlet_cell_.index_to_rank(i) != world_.rank()) {
                request.push_back(world_.isend(dirichlet_cell_.index_to_rank(i), dirichlet_cell_.index_to_binindex(i), dirichlet_cell_.cell(i)));
            }

            //++send_proc[dirichlet_cell_.index_to_rank(i)];
            assert(dirichlet_cell_.index_to_rank(i) != world_.rank());
            //if (dirichlet_cell_.index_to_rank(i) != world_.rank())
            {
                ++receive_count[dirichlet_cell_.index_to_rank(i)];
                //int a = 1;
                //MPI_Win_lock(MPI_LOCK_EXCLUSIVE, dirichlet_cell_.index_to_rank(i), 0, send_proc_win_dirichlet_);
                //MPI_Accumulate(&a, 1, MPI_INT, dirichlet_cell_.index_to_rank(i), dirichlet_cell_.index_to_rank(i), 1, MPI_INT, MPI_SUM, send_proc_win_dirichlet_);
                //MPI_Win_unlock(dirichlet_cell_.index_to_rank(i), send_proc_win_dirichlet_);
            }
            //else
                //++send_proc_dirichlet_[dirichlet_cell_.index_to_rank(i)];

            assert(tag_() != dirichlet_cell_.index_to_binindex(i));
        }
    }

    void SpCellMover::send_moving_walls(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win_wall_, int* send_proc_wall_, std::vector<int>& receive_count)
    {
        for (int i=0; i<wall_cell_.cell().size(); ++i)
        {
            for (const auto a: wall_cell_.cell(i))
            {
                assert(a.boutype() == "wall");
            }
            if (wall_cell_.index_to_rank(i) != world_.rank()) {
                request.push_back(world_.isend(wall_cell_.index_to_rank(i), wall_cell_.index_to_binindex(i), wall_cell_.cell(i)));
            }

            //++send_proc[wall_cell_.index_to_rank(i)];
            assert(wall_cell_.index_to_rank(i) != world_.rank());
            //if (wall_cell_.index_to_rank(i) != world_.rank())
            {
                ++receive_count[wall_cell_.index_to_rank(i)];

                //int a = 1;
                //MPI_Win_lock(MPI_LOCK_EXCLUSIVE, wall_cell_.index_to_rank(i), 0, send_proc_win_wall_);
                //MPI_Accumulate(&a, 1, MPI_INT, wall_cell_.index_to_rank(i), wall_cell_.index_to_rank(i), 1, MPI_INT, MPI_SUM, send_proc_win_wall_);
                //MPI_Win_unlock(wall_cell_.index_to_rank(i), send_proc_win_wall_);
            }
            //else
                //++send_proc_wall_[wall_cell_.index_to_rank(i)];

            assert(tag_() != wall_cell_.index_to_binindex(i));
        }
    }

    void SpCellMover::recv_moving_cells_self()
    {
        for (int i=0; i<movcel_.cell().size(); ++i)
        {
            if (movcel_.index_to_rank(i) == world_.rank())
            {
                arrival_cell_.insert(arrival_cell_.end(), movcel_.cell(i).begin(), movcel_.cell(i).end());
            }
        }
    }

    void SpCellMover::recv_moving_walls_self()
    {
        for (int i=0; i<wall_cell_.cell().size(); ++i)
        {
            if (wall_cell_.index_to_rank(i) == world_.rank())
            {
                arrival_wall_cell_.insert(arrival_wall_cell_.end(), wall_cell_.cell(i).begin(), wall_cell_.cell(i).end());
            }
        }
    }

    void SpCellMover::recv_moving_farfield_self()
    {
        for (int i=0; i<farfield_cell_.cell().size(); ++i)
        {
            if (farfield_cell_.index_to_rank(i) == world_.rank())
            {
                arrival_farfield_cell_.insert(arrival_farfield_cell_.end(), farfield_cell_.cell(i).begin(), farfield_cell_.cell(i).end());
            }
        }
    }

    void SpCellMover::recv_moving_empty_self()
    {
        for (int i=0; i<empty_cell_.cell().size(); ++i)
        {
            if (empty_cell_.index_to_rank(i) == world_.rank())
            {
                arrival_empty_cell_.insert(arrival_empty_cell_.end(), empty_cell_.cell(i).begin(), empty_cell_.cell(i).end());
            }
        }
    }

    void SpCellMover::recv_moving_dirichlet_self()
    {
        for (int i=0; i<dirichlet_cell_.cell().size(); ++i)
        {
            if (dirichlet_cell_.index_to_rank(i) == world_.rank())
            {
                arrival_dirichlet_cell_.insert(arrival_dirichlet_cell_.end(), dirichlet_cell_.cell(i).begin(), dirichlet_cell_.cell(i).end());
            }
        }
    }

    void SpCellMover::recv_moving_cells(const std::vector<MeshCell>& arrivals)
    {
        if (master()) return;

        for (const auto d: arrival_info_.disp())
        {
            assert(!arrivals.empty());
            if (arrival_cell_.empty()) {
                std::copy(arrivals.begin()+d.first, arrivals.begin()+d.second, std::back_inserter(arrival_cell_));
            }
            else {
                arrival_cell_.insert(arrival_cell_.end(), arrivals.begin()+d.first, arrivals.begin()+d.second);
            }

        }
    }

    void SpCellMover::reserve_arrival_cells()
    {
        if (!mergebins_)
        {
            int count = 0;
            for (int i=0; i<movcel_.cell().size(); ++i)
            {
                if (movcel_.index_to_rank(i) == world_.rank())
                {
                    ++count;
                }
            }

            arrival_cell_.reserve(arrival_info_.size() + count);
        }
        else
        {
            arrival_cell_.reserve(arrival_info_.size());
        }
    }

    void SpCellMover::recv_moving_cells(int source)
    {
        if (master()) return;

        std::vector<MeshCell> temp;
        //world_.recv(source, tag_(), temp);
        world_.recv(source, sp_->tag()(), temp);

        for (const MeshCell& mc: temp)
        {
            //if (!mc.is_ghost())
                if (mc.boutype() != "interior")
                {
                    std::cout << "is_wall: " << (mc.boutype() == "wall") << std::endl;
                    std::cout << "is_dirichlet: " << (mc.boutype() == "dirichlet") << std::endl;
                }
                assert(mc.boutype() == "interior");
                assert(sp_->aabb().do_intersect(mc.poly()));

        }
        arrival_cell_.insert(arrival_cell_.end(), temp.begin(), temp.end());
    }

    void SpCellMover::recv_wall_cells(int source, int tag)
    {
        if (master()) return;

        std::vector<MeshCell> temp;
        //std::cout << "aaaaaaaaaaaaaa - " << world_.rank() << " " << source << " " << sp_->tag()() << std::endl;
        world_.recv(source, tag, temp);
        //std::cout << "bbbbbbbbbbbbb - " << world_.rank() << std::endl;
        arrival_wall_cell_.insert(arrival_wall_cell_.end(), temp.begin(), temp.end());
    }

    void SpCellMover::recv_dirichlet_cells(int source, int tag)
    {
        if (master()) return;

        std::vector<MeshCell> temp;
        world_.recv(source, tag, temp);
        arrival_dirichlet_cell_.insert(arrival_dirichlet_cell_.end(), temp.begin(), temp.end());
    }
    
    void SpCellMover::recv_empty_cells(int source, int tag)
    {
        if (master()) return;

        std::vector<MeshCell> temp;
        world_.recv(source, tag, temp);
        arrival_empty_cell_.insert(arrival_empty_cell_.end(), temp.begin(), temp.end());
    }

    void SpCellMover::recv_farfield_cells(int source, int tag)
    {
        if (master()) return;

        std::vector<MeshCell> temp;
        world_.recv(source, tag, temp);
        arrival_farfield_cell_.insert(arrival_farfield_cell_.end(), temp.begin(), temp.end());
    }

    void ArrivalInfo::add_size(size_t s)
    {
        narrival_ += s;
    }

    void ArrivalInfo::add_disp(int begin, int end)
    {
        disp_.push_back(std::pair<int,int>(begin, end)); 
    }

    size_t ArrivalInfo::size() const
    {
        return narrival_;
    }

    const std::vector<std::pair<int,int>>& ArrivalInfo::disp() const
    {
        return disp_;
    }
    
    ArrivalInfo::ArrivalInfo(): narrival_(0) {}
}
