#ifndef COMMON_SP_H
#define	COMMON_SP_H

#include "regular_mesh.h"
//#include "movingcell.h"
#include <boost/serialization/deque.hpp>
#include <boost/serialization/vector.hpp>
#include "profiler.h"
#include "adt.h"
//#include "stencil_walk.h"

namespace Common
{
    class SpatialPartition
    {
        friend class DonorInfo;
        friend class DirectCutter;
        friend class DonorSearcher;
        friend class HoleProfiler;
        friend class SpCellMover;
        friend class Disconnecter;
        friend class ArrivalConnecter;
        friend class boost::serialization::access;

        public:

            const std::vector<RegularMesh>& rm() const;
            const Tag& tag() const;
            const AABB& aabb() const;
            const std::deque<Mesh>& mesh() const;
            const Mesh& mesh(const Tag& t) const;
            const Mesh& mesh(int i) const;
            void info() const;
            const RegularMesh& rm(const Tag& mt) const;

            void add_merge_mesh_leave_dups(const SpatialPartition& other_sp);
            void remove_dups();

            // Load
            //const std::vector<double>& mesh_load() const;
            //double load_solver(double& total_ml, double& total_nwall, double& total_nouter, double& total_npart);
            //double load(LoadEstimType load_estim_type, double& total_ml, double& total_nwall, double& total_nouter, double& total_npart);

            void reset_erase_marks();

            bool is_resident(const vec3<double>& cnt, const Outline& outline) const;
            bool is_resident(const vec3<double>& cnt) const;
            std::deque<std::vector<Point>> mesh_pts(const Outline& outline) const;
            std::deque<std::vector<Point>> mesh_pts() const;
            bool fully_resident(const MeshCell& mc) const;
            void connect_ghosts_and_neis(bool verbose=false);

            void extend_aabb(const AABB& aabb);

            void set_tag(const Tag& t);
            void set_comm(const MPI_Comm& comm);

            // Regular map
            void make_regular_maps_for_mesh();
            void make_regular_maps_for_mesh_uniproc(const std::shared_ptr<Profiler>& profiler);
            void make_regular_maps_for_mesh_serial();

            // Connect cells
            void connect_cells();
            void connect_cells_uniproc();

            // Move mesh
            void rotate_mesh(const Tag& _parent_mesh, double angle, int axis, const vec3<double>& rot_point);
            void move_mesh(const Tag& _parent_mesh, const vec3<double>& v);

            // Add mesh
            void add_mesh(const Mesh& m);
            void add_mesh(Mesh&& m);
            void add_mesh(const SpatialPartition& other_sp);

            // Merge
            void merge(SpatialPartition& other_sp);
            void merge_mesh(const SpatialPartition& other_sp);
            void merge_meshes(std::deque<Mesh>& mesh) const;
            void merge_meshes_no_check(std::deque<Mesh>& mesh) const;

            // Setter
            void set_worker_comm(const boost::mpi::communicator& comm);
            void set_profiler(std::shared_ptr<Profiler> profiler);
            void set_aabb(const AABB& aabb);

            void print_meshes() const;
            void remove_dup_cells_and_points();
            void dummy(LoadEstimType load_estim_type, double& total_ml, double& total_nwall, double& total_nouter, double& total_npart);
            void set_neighborship(const vec3<unsigned int>& global_nstripe, std::vector<unsigned int>& bin_to_proc);
            void remove_mesh(Tag mt);

            void mark_to_be_erased(const Tag& im, const Tag& ic);
            void erase_marked_cells();

            //std::deque<AABB> mesh_aabb(const Outline& outline) const;
            //std::deque<AABB> mesh_aabb() const;

            void make_cell_adt();
            void make_mesh_aabb();
            const std::deque<ADT>& cell_adt() const;
            const std::deque<AABB>& mesh_aabb() const;
            const ADT& cell_adt(int i) const;
            const AABB& mesh_aabb(int i) const;

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & mesh_;
                ar & aabb_;
                ar & tag_;
            }

        private:

            std::deque<ADT> cell_adt_;
            std::deque<AABB> mesh_aabb_;
            std::shared_ptr<Profiler> profiler_;
            Tag tag_; // same as bintag;
            boost::mpi::communicator world_;
            boost::mpi::communicator worker_comm_;
            std::vector<RegularMesh> rm_; 
            std::deque<Mesh> mesh_;
            AABB aabb_;
            //std::vector<double> mesh_load_;

            Mesh& mesh_p(const Tag& t);
            RegularMesh& rm_p(const Tag& mt);
            bool master() const;
            Tag generate_meshtag() const;
            void make_regular_map_for_mesh(const Tag& mt);
            void make_regular_map_for_mesh_uniproc(const Tag& mt, const std::shared_ptr<Profiler>& profiler);
            void make_regular_map_for_mesh_serial(const Tag& mt);
            void remove_cell(const Tag& im, const Tag& ic);
            void profstart(std::string fun);
            void profstop(std::string fun);
    };
}

#endif
