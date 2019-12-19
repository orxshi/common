#ifndef COMMON_SPSENDER_H
#define	COMMON_SPSENDER_H

#include "sp.h"
#include "profiler.h"
#include "recv_request.h"

namespace Common
{
    class SpSender
    {
        public:

            SpSender(const boost::mpi::communicator world, std::shared_ptr<Profiler> profiler);
            SpSender(const boost::mpi::communicator world);

            void send(const std::deque<Mesh>& mesh, const RegularMesh& rm, const std::vector<MeshTransferInfo>& mti, bool check_residency, std::deque<SpatialPartition>& incoming);
            void transfer_batch(const std::deque<Mesh>& mesh, const RegularMesh& rm, const std::vector<MeshTransferInfo>& mti, bool check_residency, std::deque<SpatialPartition>& incoming);
            void transfer_per_dest(const std::deque<Mesh>& mesh, const RegularMesh& rm, const std::vector<MeshTransferInfo>& mti, bool check_residency, std::deque<SpatialPartition>& incoming);

        private:

            std::shared_ptr<Profiler> profiler_;
            boost::mpi::communicator world_;
            std::vector<RecvRequest> recv_req_tup;

            void try_to_recv(std::vector<RecvRequest>& recv_req_tup, std::deque<SpatialPartition>& incoming);
            void complete_recv(const std::vector<MeshTransferInfo>& mti, std::vector<RecvRequest>& recv_req_tup, std::deque<SpatialPartition>& incoming);
            void fill_recv_req_tup(const std::vector<MeshTransferInfo>& mti, std::deque<SpatialPartition>& incoming);
            void profstart(std::string fun);
            void profstop(std::string fun);
            void prepare_sp(const std::deque<Mesh>& mesh, const RegularMesh& rm, const std::vector<MeshTransferInfo>& mti, bool check_residency, std::deque<SpatialPartition>& sps, const MeshTransferInfo& mti_, std::deque<SpatialPartition>& incoming);
            void prepare_sp_per_dest(const std::deque<Mesh>& mesh, const RegularMesh& rm, bool check_residency, SpatialPartition& sp, const MeshTransferInfo::Dest& dest);
    };
}

#endif
