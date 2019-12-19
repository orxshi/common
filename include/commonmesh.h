#ifndef COMMON_MESH_H
#define	COMMON_MESH_H

#include <vector>
#include <numeric>
#include <array>
#include <string>
#include <algorithm>
#include <boost/bimap.hpp>
#include "tag.h"
#include "geom.h"
#include "commonmeshpoint.h"
#include "commonmeshcell.h"
#include <fstream>
#include <typeinfo>        
#include "utility.h"

namespace Common
{
    typedef std::vector<MeshPoint> mpc;
    typedef std::vector<MeshCell> mcc;
    typedef std::vector<MeshFace> mfc;
    typedef std::vector<int> intvector;
    typedef boost::bimap<int, int> bimap_int;

    enum class StencilWalkResult
    {
        inside_cell = 0,
        inside_hole = 1,
        outside_mesh = 2,
    };

    class Mesh
    {    
        public: 

            Mesh() = default;
            Mesh(const Tag& tag, const Tag& parent_mesh);
            Mesh(const Tag& tag): Mesh(tag, tag) {};

            virtual void print_as_vtk(std::string file_name) const;
            void set_all_cells_as_interior();

            // face
            const mfc& face() const;

            void print_donor_info() const;
            void reset_erase_marks();
            void erase_marked_cells();
            void mark_to_be_erased(const Tag& ic);

            void bbox(vec3<double>& min_, vec3<double>& max_) const;

            // boundary
            //void add_wall_boundary(const MeshCell& mc);
            //void add_dirichlet_boundary(const MeshCell& mc);
            //void add_empty_boundary(const MeshCell& mc);
            //void add_farfield_boundary(const MeshCell& mc);
            //void add_empty_boundary(const MeshCell& mc);
            //void add_boundary(const MeshCell& mc);
            //void add_walls(std::vector<MeshCell>& cells);
            //void add_outers(std::vector<MeshCell>& cells);
            const std::vector<MeshCell>& wall_boundaries() const;
            const std::vector<MeshCell>& dirichlet_boundaries() const;
            const std::vector<MeshCell>& empty_boundaries() const;
            const std::vector<MeshCell>& farfield_boundaries() const;
            const MeshCell& empty_boundary(const Tag& t) const;
            const MeshCell& wall_boundary(const Tag& t) const;
            const MeshCell& dirichlet_boundary(const Tag& t) const;
            const MeshCell& farfield_boundary(const Tag& t) const;
            MeshCell& boundary(boundary_t btype, const Tag& t);
            const MeshCell& boundary(boundary_t btype, const Tag& t) const;
            void connect_add_bou_to_interior(Mesh& boumesh, boundary_t boutype);

            // point
            void update_points_from_cell_vertices();
            void remove_merge_duplicate_points();
            void remove_duplicate_points();
            void add_points(const std::vector<MeshCell>& cells);
            const MeshPoint& point(const Tag& ptag) const;
            const MeshPoint& point(int) const = delete;
            const std::vector<MeshPoint>& point() const;
            std::vector<Point> rawpoint() const;
            bool do_point_exist(const Tag& t) const;
            void add_point(MeshPoint p, const Tag& t, size_t size=0);
            void add_cell_only(MeshCell c, size_t size=0);
            void add_interior_nonsorted_addpoint(MeshCell mc);
            void remove_point(Tag ip);
            void remove_isolated_points();
            void sort_points();
            void shrink_points();
            void add_meshpoint(const MeshCell& mc, int rank);
            unsigned int point_index(const Tag& t) const;

            // cell
            void remove_merge_duplicate_cells();
            void remove_duplicate_cells();
            const MeshCell& cell(int) const = delete;
            const MeshCell& cell(const Tag& ctag) const;
            MeshCell& cell_p(const Tag& ctag);
            const std::vector<MeshCell>& cell() const;
            void add_interior_cell(const MeshCell& c);
            void add_interior_cell_no_check(const MeshCell& c);
            void add_interior_cell_find(const MeshCell& c);
            void add_interior_cell_nonsorted(const MeshCell& c);
            void add_interior_cell_sorted(const MeshCell& c);
            std::vector<MeshCell>::iterator add_cell_sorted(const MeshCell& c, bool& exist);
            std::vector<MeshCell>::iterator add_cell_nonsorted(const MeshCell& c);
            void add_element_nonsorted(const MeshCell& mc);
            void add_element(const MeshCell& mc);
            void add_element_no_check(const MeshCell& mc);
            void insert_cells(const std::vector<MeshCell>& cells);
            void remove_cell(const Tag& ic);
            void sort_cells();
            void connect_cells();
            void destroy_cell_hood();
            void shrink_cells();
            const MeshCell* cell_ptr(const Tag& ctag) const;

            // tag
            const Tag& tag() const;
            const Tag& parent_mesh() const;
            void set_tag(const Tag& t);
            void set_parent_mesh(const Tag& t);

            // face
            MeshFace& face_p(const Tag& ctag);
            const MeshFace& face(const Tag& ctag) const;
            const MeshFace& face(int) const = delete;
            
            // reservations
            void reserve_interior(size_t size);
            void reserve_wall(size_t size);
            void reserve_dirichlet(size_t size);
            void reserve_empty(size_t size);
            void reserve_farfield(size_t size);
            void reserve_interior_only(size_t size);

            // queries
            const MeshCell* query_sorted(const Tag& ic) const;
            const MeshPoint* query_point(const Tag& ic) const;
            const MeshCell* query(const Tag& ic) const;
            const MeshCell* query_wall(const Tag& ic) const;
            const MeshCell* query_dirichlet(const Tag& ic) const;

            // others
            //void merge_outer_to_interior(const Mesh& outer_mesh);
            void merge_batch(const Mesh& other_mesh, int rank);
            const AABB& hole_aabb() const;
            const bimap_int& dirichlet_tag_index_map() const;
            size_t npartition() const;
            void add_interiors_nonsorted_nopoint(std::vector<MeshCell>& cells);
            size_t calc_nhole_cell() const;
            void move(const vec3<double>& v);
            void rotate(double angle, int axis, const vec3<double>& rot_point);
            void merge(const Mesh& other_mesh);
            void merge_no_check(const Mesh& other_mesh);
            void remove_dup_cells_and_points();
            void simple_merge(const Mesh& other_mesh);
            void remove_parent_cells_of_vertices_of_all_cells();
            void set_parent_cell_of_vertices_of_all_cells();
            void remove_parent_cells_of_all_points();
            bool operator<(const Mesh& other) const;
            void prepare_to_remove_cell(const Tag& ic);

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & tag_;
                ar & parent_mesh_;
                ar & point_;
                ar & cell_;
                ar & wall_tag_index_map_;
                ar & dirichlet_tag_index_map_;
                ar & wall_boundaries_;
                ar & dirichlet_boundaries_;
                ar & empty_boundaries_;
                ar & farfield_boundaries_;
                ar & face_;
            }

        private:

            mcc wall_boundaries_;
            mcc dirichlet_boundaries_;
            mcc empty_boundaries_;
            mcc farfield_boundaries_;
            mcc cell_;
            mpc point_;
            mfc face_;
            Tag tag_;
            Tag parent_mesh_;
            bimap_int wall_tag_index_map_;
            bimap_int dirichlet_tag_index_map_;
            bimap_int empty_tag_index_map_;
            bimap_int farfield_tag_index_map_;
            AABB hole_aabb_;

            // point
            MeshPoint& point_p(const Tag& ptag);

            // boundary
            MeshCell& wall_boundary_p(const Tag& t);
            MeshCell& dirichlet_boundary_p(const Tag& t);
            MeshCell& empty_boundary_p(const Tag& t);
            MeshCell& farfield_boundary_p(const Tag& t);
            void remove_cell_boundary(const Tag& ic);
            void add_boundary(const MeshCell& mc, mcc& container, bimap_int& index_map);

            // cell
            void set_cell(const Tag& t, const MeshCell& c);
            void remove_cell_from_cellhood(const Tag& ic);
            void deparent_cell_from_vertices(const Tag& ic);
    };

    std::ifstream& go_to_beg_of_line(std::ifstream& file, int num);
}

#endif
