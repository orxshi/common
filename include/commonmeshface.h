#ifndef COMMON_MESHFACE_H
#define	COMMON_MESHFACE_H

#include "commonmeshpoint.h"
#include "boundary.h"

namespace Common
{
    class MeshFace
    {    
        boundary_t btype_;
        bool is_boundary_;
        std::vector<Tag> parent_cell_;
        std::vector<Tag> mesh_point_;
        Polygon face_;
        Tag tag_;
        friend class boost::serialization::access;

        public:

        MeshFace() = default;
        MeshFace(std::vector<MeshPoint>& pts);

        boundary_t btype() const;
        void set_btype(boundary_t type);
        void set_as_boundary();
        bool is_boundary() const;
        void set_tag(const Tag& t);
        const Tag& tag() const;
        void rotate(double angle, int axis, const vec3<double>& rot_point);
        void move(const vec3<double>& v);
        const Polygon& face() const;
        const std::vector<Tag>& parent_cell() const;
        const std::vector<Tag>& mesh_point() const;
        const Tag& mesh_point(int i) const;
        const Tag& parent_cell(int i) const;
        void add_parent_cell(const Tag& celltag);
        void remove_parent_cells();
        void remove_parent_cell(const Tag& ic);
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & parent_cell_;
            ar & mesh_point_;
            ar & face_;
            ar & is_boundary_;
            ar & tag_;
            ar & btype_;
        }
    };
}

#endif
