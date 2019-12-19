#ifndef COMMON_REMAP_H
#define COMMON_REMAP_H

#include "loadmap.h"
#include "spc.h"

namespace Common
{
    class Remapper
    {
        public:
            Remapper(SpatialPartitionContainer* spc, const MPI_Comm& comm, std::shared_ptr<Profiler> profiler, bool remap_batch, bool verbose);
            void remap(Loadmap& lm, bool greedy);

        private:
            bool verbose_;
            boost::mpi::communicator world_;
            std::shared_ptr<Profiler> profiler_;
            SpatialPartitionContainer* spc_;
            bool remap_batch_;

            void profstart(std::string fun);
            void profstop(std::string fun);
    };
}

#endif
