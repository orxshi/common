#ifndef COMMON_LOADMAP_H
#define	COMMON_LOADMAP_H

#include <boost/serialization/map.hpp>
#include <boost/program_options.hpp>
#include "regular_mesh.h"
#include "profiler.h"
#include "spsender.h"
#include "load_calculator.h"
#include "graph.h"

namespace Common
{
    enum class LoadPart
    {
        greedy,
        spatial
    };

    enum class RegType
    {
        aabb,
        centroid,
    };

    struct LM_Settings: Settings
    {
        void config();
    };

    class Loadmap
    {
        public:

            Loadmap(const MPI_Comm& comm, const std::shared_ptr<Profiler>& profiler, bool verbose);

            const std::map<BinRMTag, int> bintag_proc_map() const;
            const size_t nmesh() const;
            const std::vector<std::vector<BinRMTag>>& aug_bin_tag() const;
            unsigned int nworker() const;
            const RegularMesh& rm() const;
            const Mesh& mesh(const Tag& tag) const;
            //double final_dev() const;
            //void make_map(const std::vector<std::string>& file_name);
            const std::deque<Mesh>& mesh() const;
            //void merge_donor_info(std::vector<boost::mpi::request>& request);
            void move_mesh(const std::vector<vec3<double>>& v);
            void rotate_mesh(const std::vector<double>& rotation, int axis, const std::vector<vec3<double>>& rot_point);
            //void make_regular_map(const std::vector<Mesh>& mesh);
            //void read_mesh_all_procs(const std::vector<std::string>& file_name);
            //void make_map_no_read(const std::deque<Mesh>& mesh, int iteration, int nmake_map, bool adaptive, RegType reg_type);
            void make_map(int iteration, int nmake_map, bool adaptive, RegType regtype);
            void get_bin_to_proc();
            void clear();
            void distribute_mti(std::vector<MeshTransferInfo>& mti) const;
            void print(int iteration);
            size_t nbin() const;
            double refine_tol() const;
            LoadEstimType load_estim_type() const;
            void distribute_mti_2(std::vector<MeshTransferInfo>& mti) const;
            void remove_cells_from_rm();
            //void add_mesh_from_sp(const std::deque<SpatialPartition>& incoming);
            //void add_arrival_meshes(const std::deque<SpatialPartition>& incoming);
            void read_settings();
            void add_mesh(const Mesh& m);

        private:

            std::string area_rep_; // pass with ctr
            bool merge_bins_; // pass with ctr
            bool balance_;
            std::unique_ptr<Graph> graph_;
            LoadPart load_part_;
            double threshold_;
            bool verbose_;
            std::vector<Mesh> arrival_meshes;
            vec3<int> start_nstripe_;
            double refine_tol_;
            double final_dev_;
            bool uniproc_;
            std::shared_ptr<Profiler> profiler_;
            //std::vector<DonorInfo> donor_info_;
            bool printmesh_;
            bool dorefine;
            bool printlm_;
            boost::mpi::communicator world;
            unsigned int nworker_;
            int refine_limit_;
            unsigned int nrefine;
            std::unique_ptr<RegularMesh> rm_;
            std::deque<Mesh> mesh_; // use Mesh pointer instead of Mesh? We have Mesh already in Tailor.
            std::vector<std::vector<BinRMTag>> aug_bin_tag_;
            std::map<BinRMTag, int> bintag_proc_map_;
            boost::bimap<int, int> mesh_tag_index_map;
            LoadEstimType load_estim_type_;

            bool master() const;
            Tag generate_meshtag() const;
            //void add_mesh(Mesh& m);
            //void add_mesh(const std::deque<Mesh>& mesh);
            bool resolution_is_enough(BinRMTag& heaviest_bt, int nmake_map, int iteration);
            bool resolution_is_enough(BinRMTag& heaviest_bt, int nmake_map);
            //void deviation_in_load(double& dev, BinRMTag& heaviest_bt, bool& collapsed, bool verbose=false);
            void print_lm() const;
            void print_rm(const RegularMesh& rm) const;
            void print_orphan(int iteration) const;
            void print_mesh(int iteration) const;
            void refine(int nmake_map, int iteration);
            void make_regular_map();
            void make_regular_map_resident();
            void make_regular_map_uniproc();
            void make_regular_map_uniproc_resident();
            void gather_regular_maps();
            AABB max_aabb() const;
            //void remap();
            //void send_sp(const std::vector<MeshTransferInfo>& mti, std::vector<boost::mpi::request>& send_req);
            void send_sp(const std::vector<MeshTransferInfo>& mti, std::vector<RecvRequest>& recv_req_tup, std::vector<Mesh>& arrival_meshes);
            void complete_recv(const std::vector<MeshTransferInfo>& mti, std::vector<RecvRequest>& recv_req_tup, std::vector<Mesh>& arrival_meshes);
            void try_to_recv(std::vector<RecvRequest>& recv_req_tup, std::vector<Mesh>& arrival_meshes);
            //void recv_sp(const std::vector<MeshTransferInfo>& mti, std::vector<Mesh>& arrival_meshes);
            //void add_mesh(std::vector<Mesh>& mesh);
            //void add_arrival_meshes(std::vector<Mesh>& arrival_meshes);
            void reset_rm(int ns);
            void reset_rm(const vec3<int>& ns);
            void recv_arrival_meshes(const SpatialPartition& sp);
            void profstart(std::string fun);
            void profstop(std::string fun);
            void calc_load(BinRMTag& heaviest_bt, std::vector<BinRMTag>& sorted_tags, bool verbose);
            BinRMTag get_heaviest_bin(const std::vector<BinRMTag>& sorted_tags, bool verbose);
            void sort_load(std::vector<double>& sorted_load, std::vector<BinRMTag>& sorted_tags, bool verbose);
            double estimated_deviation(const std::vector<double>& aug_sorted_load);

    };

    bool greedy_partition(int nhead, const std::vector<double>& input, std::vector<double>& output, const std::vector<BinRMTag>& input_bin_tag, std::vector<std::vector<BinRMTag>>& output_bin_tag, bool verbose_);
    void fill_recv_req_tup(int rank, const std::vector<MeshTransferInfo>& mti, std::vector<RecvRequest>& recv_req_tup);
    void print_rm(const RegularMesh& rm, std::string s, const std::map<BinRMTag, int>& bintag_proc_map_);
}

#endif
