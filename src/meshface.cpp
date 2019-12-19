#include "commonmeshface.h"

namespace Common
{
    //const Tag& MeshFace::tag() const
    //{
        //return tag_;
    //}

    void MeshFace::set_as_boundary()
    {
        is_boundary_ = true;
    }

    void MeshFace::set_btype(boundary_t type)
    {
        btype_ = type;
    }
    
    bool MeshFace::is_boundary() const
    {
        return is_boundary_;
    }

    const Polygon& MeshFace::face() const
    {
        return face_;
    }

    void MeshFace::move(const vec3<double>& v)
    {
        //Point p0, p1;
        //p0.set_r(segment_.vertex(0).r(0) + v(0), segment_.vertex(0).r(1) + v(1));
        //p1.set_r(segment_.vertex(1).r(0) + v(0), segment_.vertex(1).r(1) + v(1));
        //segment_ = Segment(p0, p1);
        face_.move_points(v);
    }

    void MeshFace::rotate(double angle, int axis, const vec3<double>& rot_point)
    {
        face_.rotate_points(angle, axis, rot_point);
    }

    const Tag& MeshFace::mesh_point(int i) const
    {
        assert(i >= 0);
        assert(i < mesh_point_.size());

        return mesh_point_[i];
    }

    const std::vector<Tag>& MeshFace::mesh_point() const
    {
        return mesh_point_;
    }

    void MeshFace::set_tag(const Tag& t)
    {
        tag_ = t;
    }

    boundary_t MeshFace::btype() const
    {
        return btype_;
    }

    MeshFace::MeshFace(std::vector<MeshPoint>& pts): is_boundary_(false)
    //MeshFace::MeshFace(std::vector<MeshPoint>& pts)
    {
        assert(!pts.empty());
        std::vector<Point> rawpts;
        rawpts.reserve(pts.size());

        for (const MeshPoint& mp: pts)
        {
            mesh_point_.push_back(mp.tag());
            rawpts.push_back(mp.p());
        }

        assert(!rawpts.empty());
        face_ = Polygon(rawpts);
        assert(!face_.vertex().empty());
    }

    void MeshFace::remove_parent_cell(const Tag& ic)
    {
        assert(ic.isvalid());
        parent_cell_.erase(std::remove(parent_cell_.begin(), parent_cell_.end(), ic), parent_cell_.end());
    }

    /*const Tag& MeshFace::parent_mesh() const
    {
        return parent_mesh_;
    }*/

    void MeshFace::remove_parent_cells()
    {
        parent_cell_.clear();
    }

    /*void MeshFace::set_tag(const Tag& t)
    {
        tag_ = t;
    }
    void MeshFace::set_parent_mesh(const Tag& ptag)
    {
        parent_mesh_ = ptag;
    }*/

    const std::vector<Tag>& MeshFace::parent_cell() const
    {
        return parent_cell_;
    }

    const Tag& MeshFace::parent_cell(int i) const
    {
        return parent_cell_[i];
    }

    const Tag& MeshFace::tag() const
    {
        return tag_;
    }

    void MeshFace::add_parent_cell(const Tag& celltag)
    {
        assert(celltag.isvalid());
        parent_cell_.push_back(celltag);
    }
}
