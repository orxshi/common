#ifndef COMMON_SPC_H
#define	COMMON_SPC_H

#include <boost/mpi.hpp>
#include "spsender.h"
#include "commonsettings.h"

namespace Common
{
    struct SPC_Settings: Settings
    {
        void config();
    };

    class SpatialPartitionContainer
    {
        friend class DonorInfo;
        friend class DonorSearcher;
        friend class HoleProfiler;
        friend class CellExchanger;
        friend class Disconnecter;
        friend class ArrivalConnecter;

        public:
            //void make_outline(std::vector<std::vector<vec3<double>>> allpts);
            const std::map<BinRMTag, int>& bintag_proc_map() const;

            SpatialPartitionContainer(const MPI_Comm& comm, std::shared_ptr<Profiler> profiler, bool verbose, const std::map<BinRMTag, int>& bintag_proc_map, int nmesh);
            SpatialPartitionContainer(const MPI_Comm& comm, std::shared_ptr<Profiler> profiler, bool verbose, LoadEstimType let, int nmesh);
            ~SpatialPartitionContainer();

            //void output_load(int iteration);
            
            //void reduce_sp(std::deque<SpatialPartition>& incoming);

            void make_mesh_aabbs();

            // Connect
            void connect_cells();
            void connect_cells_uniproc();

            // Getter
            const std::vector<SpatialPartition>& sp() const;
            void info() const;
            const RegularMesh& global_rm() const;

            // Make rm
            void make_regular_maps();
            void make_regular_maps_uniproc();

            // Move mesh
            void rotate_meshblocks(const Tag& _parent_mesh, double ang, int axis, const vec3<double>& rot_axis);
            void move_meshblocks(const Tag& _parent_mesh, const vec3<double>& v);

            // Remap
            void push_sp(std::deque<SpatialPartition>& sp, const std::vector<bool>& cond);
            void add_sp(SpatialPartition& other_sp);
            void add_sp(const std::deque<SpatialPartition>& sp);
            void add_sp(const std::deque<SpatialPartition>& sp, std::vector<bool>& cond);
            void remap_uniproc(const std::deque<Mesh>& mesh);

            void merge_sp(const std::deque<SpatialPartition>& sp);
            void remove_dup_cells_and_points();
            double dev() const;
            void merge_sp_meshes(std::deque<Mesh>& mesh);
            void print_all_meshes_in_partitions();
            void print_all_meshes_in_partitions_uniproc();
            bool load();
            void set_global_rm(const RegularMesh& rm);

            //void reduce_sp();

            const std::vector<Outline>& outline() const;

        private:
            std::vector<Outline> outline_;
            bool verbose_;
            bool remap_batch_;
            RegularMesh global_rm_; // global rm / load map rm / bins are empty of cells / only used as a map.
            //double load_;
            bool print_mesh_;
            double dev_;
            std::shared_ptr<Profiler> profiler_;
            boost::mpi::communicator world_;
            boost::mpi::communicator worker_comm;
            std::vector<SpatialPartition> sp_;
            std::map<BinRMTag, int> bintag_proc_map_;
            int nmesh_;
            //std::vector<SpCellMover> sp_cell_mover_;

            bool master() const;
            void read_settings();
            void profstart(std::string fun);
            void profstop(std::string fun);
    };

    bool duplicate_exist(const SpatialPartitionContainer& spc);
}

#endif
