#ifndef COMMON_MESHCELL_H
#define COMMON_MESHCELL_H

#include "commonmeshface.h"
#include <boost/serialization/utility.hpp>
#include "array.h"

#define NPOINT 8

namespace Common
{
    class MeshCell
    {
        public:

            using bouarray = Array<Tag, 6>;
            using facearray = Array<MeshFace, 6>;
            using neiarray = Array<Tag, 6>;
            using pointarray = Array<MeshPoint, NPOINT>;

            MeshCell(const Tag& tag, const Tag& parent_mesh): tag_(tag), parent_mesh_(parent_mesh), OGA_cell_type_(OGA_cell_type_t::undefined) {}
            MeshCell(): MeshCell(Tag(-1), Tag(-1)) {}
            MeshCell(const Tag& tag, const Tag& parent_mesh, const std::vector<MeshPoint>& point, boundary_t btype=boundary_t("undefined"), Shape shape=Shape::undef);
            MeshCell(const std::vector<MeshPoint>& point, boundary_t btype=boundary_t("undefined")): MeshCell(Tag(-1), Tag(-1), point, btype) {}

            boundary_t boutype() const;

            void set_oga_cell_type(OGA_cell_type_t t);

            void unmark_to_be_erased();
            void mark_to_be_erased();
            bool erase() const;

            // neighborhood
            void remove_self_from_neighbors(const Tag& t);
            void add_pnei(const Tag& celltag);
            void remove_neighbor(const Tag& t);
            void remove_all_neighbors();
            const neiarray& pnei() const;
            const Tag& pnei(int i) const;

            // query
            //bool is_wall() const;
            //bool is_dirichlet() const;
            //bool is_interior() const;
            bool near_boundary() const;
            //bool near_wall() const;
            //bool near_dirichlet() const;
            //const bool is_partition() const;

            // face
            void make_faces(Shape shape);
            void deparent_self_from_faces();
            std::vector<std::vector<Tag>> face_vertex();
            const facearray& face() const;
            facearray& face_p();
            void deparent_cell_from_faces(const Tag& tag);

            // vertex/point
            void add_vertex(const MeshPoint& p, const Tag& pointtag);
            void deparent_self_from_vertices();
            void remove_parent_cells_of_vertices();
            void set_parent_cell_of_vertices();
            const pointarray& point() const;
            const MeshPoint& point(int i) const;
            std::vector<Point> geom_point() const;
            void rotate_points(double ang, int axis, const vec3<double>& rot_axis);
            void move_points(const vec3<double>& final_loc);
            size_t npoint() const;

            // boundary
            void set_interior_boundary(const Tag& t);
            void add_boundary(const Tag& t, boundary_t boutype);
            const bouarray& wall_boundary() const;
            const bouarray& dirichlet_boundary() const;
            const bouarray& empty_boundary() const;
            const bouarray& farfield_boundary() const;
            const Tag& interior_boundary() const;

            // tag
            void set_tag(const Tag& t);
            void set_parent_mesh(const Tag& celltag);
            const Tag& parent_mesh() const;
            const Tag& tag() const;
            void set_point_tag(int i, const Tag& t);

            // other
            const Polyhedron& poly() const;
            const std::pair<Tag, Tag>& donor() const;
            OGA_cell_type_t oga_cell_type() const;
            //void set_partition(bool p);
            void set_boutype(boundary_t btype);
            void merge(const MeshCell& other);
            bool operator<(const MeshCell& other) const;
            bool operator==(const MeshCell& other) const;

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & tag_;
                ar & parent_mesh_;
                ar & donor_;
                ar & pnei_;
                ar & point_;
                ar & face_;
                ar & poly_;
                ar & OGA_cell_type_;
                ar & wall_boundary_;
                ar & dirichlet_boundary_;
                ar & dirichlet_boundary_;
                ar & empty_boundary_;
                ar & farfield_boundary_;
                ar & boutype_;
                //ar & is_wall_;
                //ar & is_dirichlet_;
                //ar & is_interior_;
            }

        private:
            boundary_t boutype_;
            bool erase_;
            bouarray wall_boundary_;
            bouarray dirichlet_boundary_;
            bouarray empty_boundary_;
            bouarray farfield_boundary_;
            Tag interior_boundary_; // only for boundaries.
            Tag tag_;
            Tag parent_mesh_;
            std::pair<Tag, Tag> donor_;
            neiarray pnei_; // principal neighbors. // dynamic
            pointarray point_;
            facearray face_; //dynamic
            Polyhedron poly_;
            OGA_cell_type_t OGA_cell_type_;

            friend class boost::serialization::access;
            friend class Mesh;
    };
}

#endif
