#include "spsender.h"

namespace Common
{
    void SpSender::profstart(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->start(fun);
    }

    void SpSender::profstop(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->stop(fun);
    }

    SpSender::SpSender(const boost::mpi::communicator world, std::shared_ptr<Profiler> profiler): world_(world, boost::mpi::comm_attach), profiler_(profiler)
    {
    }

    SpSender::SpSender(const boost::mpi::communicator world): world_(world, boost::mpi::comm_attach)
    {
        profiler_ = nullptr;
    }

    void SpSender::fill_recv_req_tup(const std::vector<MeshTransferInfo>& mti, std::deque<SpatialPartition>& incoming)
    {
        for (const MeshTransferInfo& mti_: mti)
        {
            for (const MeshTransferInfo::Dest& dest: mti_.dest())
            {
                if (world_.rank() == dest.rank_)
                {
                    recv_req_tup.push_back(RecvRequest(mti_.source(), dest.sp_tag_));
                }
            }
        }

        //incoming.reserve(recv_req_tup.size());
    }

    void SpSender::try_to_recv(std::vector<RecvRequest>& recv_req_tup, std::deque<SpatialPartition>& incoming)
    {
        auto status = world_.iprobe();
        if (status)
        {
            auto it = std::find_if(recv_req_tup.begin(), recv_req_tup.end(), [&](const RecvRequest& r)
                    {
                    return (r.source == status->source() && pairing_function(r.sp_tag.rmtag()(), r.sp_tag.bintag()()) == status->tag());
                    });
            assert(it != recv_req_tup.end());
            SpatialPartition sp;
            world_.recv(status->source(), status->tag(), sp);
            incoming.push_back(sp); // uncomment
            //fnc_(sp);
            it->received = true;
        }
    }

    void SpSender::complete_recv(const std::vector<MeshTransferInfo>& mti, std::vector<RecvRequest>& recv_req_tup, std::deque<SpatialPartition>& incoming)
    {
        //if (master()) return;
        if (world_.rank() == 0) return;

        for (const MeshTransferInfo& mti_: mti)
        {
            for (const MeshTransferInfo::Dest& dest: mti_.dest())
            {
                if (world_.rank() == dest.rank_)
                {
                    SpatialPartition sp;
                    auto it = std::find_if(recv_req_tup.begin(), recv_req_tup.end(), [&](const RecvRequest& r)
                            {
                            return (r.source == mti_.source() && r.sp_tag == dest.sp_tag_);
                            });
                    assert(it != recv_req_tup.end());
                    if (!it->received)
                    {
                        world_.recv(mti_.source(), pairing_function(it->sp_tag.rmtag()(), it->sp_tag.bintag()()), sp);
                        incoming.push_back(std::move(sp)); // uncommnet
                        //fnc_(sp);
                        it->received = true;
                    }
                }
            }
        }
    }

    void SpSender::prepare_sp(const std::deque<Mesh>& mesh, const RegularMesh& rm, const std::vector<MeshTransferInfo>& mti, bool check_residency, std::deque<SpatialPartition>& sps, const MeshTransferInfo& mti_, std::deque<SpatialPartition>& incoming)
    {
        fill_recv_req_tup(mti, incoming);

        for (const MeshTransferInfo::Dest& dest: mti_.dest())
        {
            SpatialPartition sp;
            std::deque<Mesh> mb = dest.make_meshes(mesh, rm.bin(dest.sp_tag_).aabb(), check_residency);

            const Bin& b = rm.bin(dest.sp_tag_);

            sp.set_aabb(b.aabb());
            sp.set_tag(b.tag()); 
            assert(sp.tag().isvalid());
            assert(b.aabb().min(0) != b.aabb().max(0));
            for (Mesh& m: mb)
            {
                assert(m.parent_mesh().isvalid());
                if (!m.cell().empty())
                {
                    sp.add_mesh(m);
                }
            }

            assert(!sp.mesh().empty());

            sps.push_back(std::move(sp));
            assert(!sps.back().aabb().faces().empty());
        }
    }

    void SpSender::prepare_sp_per_dest(const std::deque<Mesh>& mesh, const RegularMesh& rm, bool check_residency, SpatialPartition& sp, const MeshTransferInfo::Dest& dest)
    {
        std::deque<Mesh> mb = dest.make_meshes(mesh, rm.bin(dest.sp_tag_).aabb(), check_residency);

        const Bin& b = rm.bin(dest.sp_tag_);

        sp.set_aabb(b.aabb());
        assert(!sp.aabb().faces().empty());
        sp.set_tag(b.tag()); 
        assert(sp.tag().isvalid());
        assert(b.aabb().min(0) != b.aabb().max(0));
        int ncell = 0;
        for (Mesh& m: mb)
        {
            ncell += m.cell().size();
            assert(m.parent_mesh().isvalid());
            if (!m.cell().empty())
            {
                sp.add_mesh(m);
            }
        }

        assert(!sp.mesh().empty());
    }

    void SpSender::transfer_per_dest(const std::deque<Mesh>& mesh, const RegularMesh& rm, const std::vector<MeshTransferInfo>& mti, bool check_residency, std::deque<SpatialPartition>& incoming)
    {
        std::vector<MeshTransferInfo>::const_iterator mti_;

        if (world_.rank() != 0)
        {
            int count = std::count_if(mti.begin(), mti.end(), [&](const MeshTransferInfo& tmti){return tmti.source() == world_.rank();});
            assert(count <= 1);
            mti_ = std::find_if(mti.begin(), mti.end(), [&](const MeshTransferInfo& tmti){return tmti.source() == world_.rank();});
            if (mti_ == mti.end()) {
                return;
            }
        }

        fill_recv_req_tup(mti, incoming);

        if (world_.rank() != 0)
        {
            for (const MeshTransferInfo::Dest& dest: mti_->dest())
            {
                SpatialPartition sp;
                prepare_sp_per_dest(mesh, rm, check_residency, sp, dest);
                boost::mpi::request send_req = world_.isend(dest.rank_, pairing_function(dest.sp_tag_.rmtag()(), dest.sp_tag_.bintag()()), sp);
                while (true)
                {
                    if (send_req.test()) {
                        break;
                    }

                    try_to_recv(recv_req_tup, incoming);
                }
                //send_req.wait(); //why? you already did test().
            }

            complete_recv(mti, recv_req_tup, incoming);
        }
    }

    void SpSender::transfer_batch(const std::deque<Mesh>& mesh, const RegularMesh& rm, const std::vector<MeshTransferInfo>& mti, bool check_residency, std::deque<SpatialPartition>& incoming)
    {
        std::vector<boost::mpi::request> send_req;
        std::deque<SpatialPartition> sps;
        std::vector<MeshTransferInfo>::const_iterator mti_;
        //if (world_.rank() == 0) return;

        //profstart("spsender-count");
        if (world_.rank() != 0)
        {
            int count = std::count_if(mti.begin(), mti.end(), [&](const MeshTransferInfo& tmti){return tmti.source() == world_.rank();});
            assert(count <= 1);
            mti_ = std::find_if(mti.begin(), mti.end(), [&](const MeshTransferInfo& tmti){return tmti.source() == world_.rank();});
            if (mti_ == mti.end()) {
                return;
            }
        }
        //profstop("spsender-count");

        //profstart("spsender-prepare");
        if (world_.rank() != 0)
        {
            prepare_sp(mesh, rm, mti, check_residency, sps, *mti_, incoming);
            assert(!sps.empty());
        }
        //profstop("spsender-prepare");

        //profstart("spsender-send");
        if (world_.rank() != 0)
        {
            for (int i=0; i<mti_->dest().size(); ++i)
            {
                const MeshTransferInfo::Dest& dest = mti_->dest()[i];

                send_req.push_back(world_.isend(dest.rank_, pairing_function(dest.sp_tag_.rmtag()(), dest.sp_tag_.bintag()()), sps[i]));
            }
        }
        //profstop("spsender-send");

        //profstart("spsender-complete");
        if (world_.rank() != 0)
        {
            complete_recv(mti, recv_req_tup, incoming);
        }
        //profstop("spsender-complete");

        //profstart("spsender-wait");
        if (world_.rank() != 0)
        {
            boost::mpi::wait_all(send_req.begin(), send_req.end());
        }
        //profstop("spsender-wait");
    }

    void SpSender::send(const std::deque<Mesh>& mesh, const RegularMesh& rm, const std::vector<MeshTransferInfo>& mti, bool check_residency, std::deque<SpatialPartition>& incoming)
    {
        if (world_.rank() == 0) return;
        if (world_.rank() != 0)
        {
            //profstart("fill-recv");
            fill_recv_req_tup(mti, incoming);
            //profstop("fill-recv");

            //profstart("spsender-1");
            int count = std::count_if(mti.begin(), mti.end(), [&](const MeshTransferInfo& tmti){return tmti.source() == world_.rank();});
            assert(count <= 1);
            auto mti_ = std::find_if(mti.begin(), mti.end(), [&](const MeshTransferInfo& tmti){return tmti.source() == world_.rank();});
            //profstop("spsender-1");
            if (mti_ == mti.end()) {
                return;
            }

            for (const MeshTransferInfo::Dest& dest: mti_->dest())
            {
                SpatialPartition sp;
                //profstart("dest-make-meshes");
                std::deque<Mesh> mb = dest.make_meshes(mesh, rm.bin(dest.sp_tag_).aabb(), check_residency);
                //profstop("dest-make-meshes");

                //profstart("spsender-isend");
                const Bin& b = rm.bin(dest.sp_tag_);

                sp.set_aabb(b.aabb());
                sp.set_tag(b.tag()); 
                //sp.set_global_rm(rm); 
                assert(sp.tag().isvalid());
                assert(b.aabb().min(0) != b.aabb().max(0));
                for (Mesh& m: mb)
                {
                    assert(m.parent_mesh().isvalid());
                    if (!m.cell().empty())
                    {
                        sp.add_mesh(m);
                    }
                }

                assert(!sp.mesh().empty());

                boost::mpi::request req = world_.isend(dest.rank_, pairing_function(dest.sp_tag_.rmtag()(), dest.sp_tag_.bintag()()), sp);
                while (true)
                {
                    if (req.test()) {
                        break;
                    }

                    try_to_recv(recv_req_tup, incoming);
                }
                //profstop("spsender-isend");
            }

            //profstart("spsender-complete");
            complete_recv(mti, recv_req_tup, incoming);
            //profstop("spsender-complete");
        }
    }
}
