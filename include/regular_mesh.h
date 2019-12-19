#ifndef COMMON_REGULAR_MESH_H
#define COMMON_REGULAR_MESH_H

#include <map>
#include "commoncore.h"
#include "profiler.h"
#include "rm_bin.h"
#include <fstream>
#include <numeric>

#define RM_NNEI 9

namespace Common
{
    enum class RelativePosition
    {
        undefined = -1,
        sw = 0,
        sc = 1,
        se = 2,
        cw = 3,
        cc = 8,
        ce = 4,
        nw = 5,
        nc = 6,
        ne = 7,
    };

    bool in_range(const RegularMeshIndex& index, const RegularMeshIndex& global_min, const RegularMeshIndex& global_max);
    unsigned int mat_to_vec_index(const RegularMeshIndex& index, const vec3<unsigned int>& nstripe);
    RelativePosition get_relative_position(int index, int targetindex, const vec3<unsigned int>& nstripe, size_t maxbinsize, bool& exact);
    int get_binindex(unsigned int index, RelativePosition rp, const vec3<unsigned int>& nstripe, unsigned int binsize);
    vec3<double> llcoor(const RegularMeshIndex& index, const vec3<double>& aabb_min, const vec3<double>& h);

    class RegularMesh
    {
        public:

            RegularMesh();
            RegularMesh(const Tag& rmtag);
            RegularMesh(const RegularMesh& other);
            RegularMesh& operator=(const RegularMesh& other);

            void flood(int start, const AABB& aabb);

            void move(const vec3<double>& v);
            void rotate(double ang, int axis, const vec3<double>& rot_point);

            // Getters
            const Tag& tag() const;
            size_t size() const;
            const std::map<Tag, RegularMesh*>& rmtag_address_map() const;
            const Bin& bin(int row, int col, int depth) const;
            const Bin& bin(const RegularMeshIndex& ind) const;
            const Bin& bin(int i) const;
            const Bin& bin(const Tag& t) const;
            const Bin& bin(const BinRMTag& bt) const;
            RegularMeshIndex global_index(const RegularMeshIndex& i) const;
            const RegularMeshIndex& global_min() const;
            const RegularMeshIndex& global_max() const;
            const vec3<int>& nstripe() const;
            int nstripe(int i) const;
            const std::vector<Bin>& bin() const;
            const vec3<double>& h() const;
            double h(int i) const;
            const AABB& aabb() const;

            // refine
            void refine(const std::vector<Mesh>& meshes, int rank);
            void refine_adaptive(const std::deque<Mesh>& meshes, const BinRMTag& heaviest_bt, int rank);
            void refine_adaptive(const std::deque<Mesh>& meshes, const BinRMTag& heaviest_bt, int rank, bool uniproc, int& new_rmtag, int& new_bt);

            // register
            void register_cell(const MeshCell& cell, const Tag& cell_tag, const Tag& mesh_tag, bool calcload, int rank, bool adaptive, double& dur1, double& dur2, double& dur3, double& dur4, double& dur5);
            void register_resident_cell(const MeshCell& cell, const Tag& cell_tag, const Tag& mesh_tag, bool calcload, int rank);
            void register_mesh(const Mesh& mesh, bool calcload, int rank, bool adaptive, const std::shared_ptr<Profiler>& profiler);
            void register_overlapping_mesh(const Mesh& mesh, bool calcload, int rank, bool adaptive);
            void register_resident_mesh(const Mesh& mesh, bool calcload, int rank);
            void register_bincells(const std::vector<BinCell>& bincell, const std::deque<Mesh>& meshes);

            BinRMTag index_to_bin(size_t i) const;
            //void gather_load(const boost::mpi::communicator& comm, std::vector<double>& load_global, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
            //void gather_load(std::vector<double>& load_global, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
            void calc_aabb_max();
            bool extend_aabb(const AABB& other_aabb);
            void insert_to_rmtag_address_map(const Tag& t, RegularMesh* rmp);
            void set_rmtag_address_map(const std::map<Tag, RegularMesh*>& map);
            void update_address(std::map<Tag, RegularMesh*>& map);
            void update_address();
            void set_tag(const Tag& t);
            //std::vector<BinRMTag> sort_bin_tags_based_on_load(const boost::mpi::communicator& comm, std::vector<double>& load_global, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
            void get_bin_tags(std::vector<BinRMTag>& tag) const;
            void calc_h();
            void clear_bincells();
            int get_binindex(unsigned int index, RelativePosition rp);
            std::vector<int> pneis_of_bin(unsigned int index);
            void set_bincell(const RegularMeshIndex& ind, const std::vector<BinCell>& bc);
            void set_global_min(const RegularMeshIndex& i);
            void set_global_max(const RegularMeshIndex& i);
            void set_global_min(int i, int j, int k);
            void set_global_max(int i, int j, int k);
            void set_nstripe(int x, int y, int z);
            void set_nstripe(const vec3<int>& ns);
            void set_h(const vec3<double>& v);
            void set_h(double v0, double v1, double v2);
            void set_aabb(const AABB& aabb);
            //void set_aabb_min(const vec3<double>& m);
            //void set_aabb_max(const vec3<double>& m);
            vec3<double> centroid(const RegularMeshIndex&) const;
            vec3<double> llcoor(const RegularMeshIndex& index) const;
            bool get_index_regular(const vec3<double>& cv, std::vector<RegularMeshIndex>& index, bool closest) const;
            bool get_bintag_adaptive(const vec3<double>& cv, std::vector<BinRMTag>& tag, bool closest) const;
            void get_bintag_adaptive(const AABB& aabb, std::vector<BinRMTag>& tag, int& nlevel) const;
            void get_bintag_adaptive_(const AABB& aabb, std::vector<BinRMTag>& tag, int& nlevel) const;
            void print(std::string file_name, const std::map<BinRMTag, int>& bintag_proc_map) const;
            void print(std::string file_name) const;
            void merge(const RegularMesh& other);
            bool point_inside_bin(const RegularMeshIndex& ind, const vec3<double>& p);
            void insert_bins(int mintag, int mesh_load_size);
            void set_props(const Mesh& mesh, const AABB& forced_min_aabb, int _nstripe=0);
            void info() const;
            bool validate(const RegularMeshIndex& global_ind, RegularMeshIndex& local_ind) const;
            void calc_step_length();
            void clear_cells();
            bool operator<(const RegularMesh& other) const;
            bool operator==(const RegularMesh& r) const;

        private:

            Tag tag_;
            vec3<int> nstripe_;
            vec3<double> h_;
            AABB aabb_;
            std::vector<Bin> bin_;
            RegularMeshIndex global_min_;
            RegularMeshIndex global_max_;
            std::map<Tag, RegularMesh*> rmtag_address_map_;
            friend class boost::serialization::access;
            friend class HoleMap;

            Bin& bin_p(int row, int col, int depth);
            Bin& bin_p(const RegularMeshIndex& ind);
            Bin& bin_p(const Tag& t);
            Bin& bin_p(const BinRMTag& bt);
            bool is_resident(const vec3<double>& cnt) const;
            int get_new_bt_p() const;
            int get_new_bt() const;
            int get_new_rmtag() const;
            //void gather_load(std::vector<double>& load_local, int& j, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
            void size(size_t& s) const;
            BinRMTag index_to_bin(size_t& s, size_t i) const;

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & nstripe_;
                ar & bin_;
                ar & tag_;
                ar & rmtag_address_map_;
                ar & aabb_;
                ar & h_;
            }
    };

    bool get_index_regular(const vec3<double>& cv, std::vector<RegularMeshIndex>& index, const vec3<double>& aabb_min, const vec3<double>& h, const vec3<int>& nstripe, bool closest);
    Tag closest_non_empty_bin(const RegularMeshIndex& ind, const RegularMesh& rm, const Mesh& m);
}

#endif
