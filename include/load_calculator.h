#ifndef COMMON_LOAD_CALCULATOR_H
#define COMMON_LOAD_CALCULATOR_H

#include "spc.h"
#include "rm_bin.h"

namespace Common
{
    class LoadCalculator
    {
        public:

            LoadCalculator(LoadEstimType type, std::string area_rep_, bool merge_bins_, const MPI_Comm& comm);
            bool load(const SpatialPartitionContainer& spc, double refine_tol, int iteration) const;
            std::vector<BinRMTag> sorted_bin_tags(const RegularMesh& rm, std::vector<double>& load_global, std::vector<double>& load_global_r, const std::deque<Mesh>& mesh);
            void gather_load(const RegularMesh& rm, std::vector<double>& load_global, std::vector<double>& load_global_r, const std::deque<Mesh>& mesh);
            double load(const Bin& bin, const std::deque<Mesh>& mesh) const;
            double load_r(const Bin& bin, const std::deque<Mesh>& mesh) const;

        private:

            LoadEstimType type_;
            boost::mpi::communicator world_;
            //std::shared_ptr<Settings> settings_;
            std::string area_rep_;
            bool merge_bins_;

            //double area_load_alpha(const Bin& bin, const std::deque<Mesh>& mesh) const;
            //double area_load_convex(const Bin& bin, const std::deque<Mesh>& mesh) const;
            double area_load(const Bin& bin, const std::deque<Mesh>& mesh) const;
            //double area_load(const SpatialPartition& sp, const Outline& outline) const;
            double area_load(const SpatialPartition& sp) const;
            //double area_load_alpha(const SpatialPartition& sp, const Outline& outline) const;
            //double area_load_alpha(const SpatialPartition& sp) const;
            //double area_load_convex(const SpatialPartition& sp, const Outline& outline) const;
            //double area_load_convex(const SpatialPartition& sp) const;
            double hybrid_load(const Bin& bin) const;
            //double hybrid_load(const SpatialPartition& sp, const Outline& outline) const;
            double hybrid_load(const SpatialPartition& sp) const;
            double solver_load(const Bin& bin) const;
            //double solver_load(const SpatialPartition& sp, const Outline& outline) const;
            double solver_load(const SpatialPartition& sp) const;
            double minmesh_load(const Bin& bin) const;
            //double minmesh_load(const SpatialPartition& sp, const Outline& outline) const;
            //double load(const SpatialPartition& sp, const Outline& outline) const;
            double load(const SpatialPartition& sp) const;
            bool master() const;
            void gather_load(const RegularMesh& rm, std::vector<double>& load_local, std::vector<double>& load_local_r, int& j, const std::deque<Mesh>& mesh);
    };
}

#endif
