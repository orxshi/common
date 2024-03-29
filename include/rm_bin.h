#ifndef COMMON_REGULAR_MESH_BIN_H
#define COMMON_REGULAR_MESH_BIN_H

#include "tag.h"
#include "rm_index.h"
#include "commonmesh.h"
#include "mesh_transfer_info.h"

namespace Common
{
    class RegularMesh;

    class BinCell
    {
        public:

            const Tag& mesh() const;
            const Tag& cell() const;
            int proc() const;
            //const Tag& tag() const;
            //void set_tag(const Tag& tag);
            void set_mesh(const Tag& tag);
            void set_cell(const Tag& tag);
            void set_proc(int proc);
            bool operator<(const BinCell& other) const;

        private:

            //Tag tag_;
            Tag cell_;
            Tag mesh_;
            int proc_;
            friend class boost::serialization::access;

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & cell_;
                ar & mesh_;
                ar & proc_;
            }
    };

    class Bin
    {
        public:

            Bin(const Tag& t, const Tag& parent_rm, const RegularMeshIndex& index, const AABB& aabb, int mesh_load_size=0);
            Bin(const Tag& t);
            //Bin(): load_(0.) {}
            Bin(): rm_(nullptr) {}
            ~Bin();
            Bin(const Bin& other);
            Bin& operator=(const Bin& other);

            void set_load(double l);
            double load() const;

            void move(const vec3<double>& v);
            void rotate(double ang, int axis, const vec3<double>& rot_point);

            std::deque<std::vector<Point>> mesh_pts(const std::deque<Mesh>& mesh) const;
            const boost::bimap<int, int>& mesh_tag_index_map() const;
            void register_cells_to_rm(const std::deque<Mesh>& meshes);
            void init_rm(const Bin& other);
            //void init_rm(int mintag);
            void init_rm(int mintag, const Tag& rmtag, const vec3<int>& stripe);
            const Tag& tag() const;
            void clear_cells();
            //boost::shared_ptr<RegularMesh> rm() const;
            RegularMesh* rm() const;
            void merge(const Bin& other);
            const boost::bimap<int, int>& mesh_list() const;
            bool is_resident(const vec3<double>& cnt) const;
            void set_aabb(const AABB& aabb);
            //std::deque<Mesh> group_mesh(const std::deque<Mesh>& mesh) const;
            void group_mesh(const std::deque<Mesh>& mesh, std::deque<Mesh>& mb) const;
            //void append_to_cell(const std::vector<BinCell>& other);
            //bool query_bincell(const Tag& im, const Tag& ic) const;
            void resize_mesh_load(size_t size);
            const AABB& aabb() const;
            const std::vector<double>& mesh_load() const;
            void set_index(const RegularMeshIndex& ind);
            const RegularMeshIndex& index() const;
            const std::vector<BinCell>& cell() const;
            const BinCell& cell(int i) const;
            void set_cell(const std::vector<BinCell>& c);
            double load(const std::deque<Mesh>& mesh, LoadEstimType load_estim_type);
            double load_with_area(const std::deque<Mesh>& mesh);
            double load_solver(const std::deque<Mesh>& mesh);
            void copy_bincell(const BinCell& c);
            void add_bincell(const Tag& ctag, const Tag& mtag, int proc);
            void increment_mesh_load(const Tag& t, int load=-1);
            int mesh_load(const Tag& t) const;
            void prepare_transfer_info(int dest, std::vector<MeshTransferInfo>& mti, std::vector<int>& mti_proc, std::vector<int>& dest_rank, int nbin) const;
            //void prepare_transfer_info_with_mesh(int dest, std::vector<MeshTransferInfo>& mti, std::vector<int>& mti_proc, std::vector<int>& dest_rank, const std::deque<Mesh>& mesh) const;
            const Tag& parent_rm() const;
            double load_without_calc() const;
            void reserve(size_t size);
            void shrink();
            bool operator<(const Bin& other) const;
            bool operator==(const Bin& other) const;
            std::deque<AABB> mesh_aabb(const std::deque<Mesh>& mesh) const;

        private:

            Tag tag_;
            Tag parent_rm_;
            RegularMeshIndex rmindex_;
            double load_;
            std::vector<double> mesh_load_; /*!< Number of cells belonging to meshes.*/
            //MeshLoad mesh_load_;
            std::vector<BinCell> cell_;
            boost::bimap<int, int> mesh_tag_index_map_;
            //boost::bimap<BinCell, int> bincell_index_map;
            //std::map<BinCell, int> bincell_index_map;
            AABB aabb_;
            RegularMesh* rm_;
            friend class boost::serialization::access;
            friend class HoleMap;

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & rmindex_;
                ar & load_;
                ar & cell_;
                ar & mesh_load_;
                ar & mesh_tag_index_map_;
                ar & rm_;
                ar & tag_;
                ar & parent_rm_;
                ar & aabb_;
            }
    };
}

#endif
