#include <commonmeshcell.h>

namespace Common
{
    bool MeshCell::operator<(const MeshCell& other) const
    {
        return tag_ < other.tag();
    }

    bool MeshCell::operator==(const MeshCell& other) const
    {
        return tag_ == other.tag();
    }

    void MeshCell::merge(const MeshCell& other)
    {
    }

    size_t MeshCell::npoint() const
    {
        return point_.size();
    }

    const Tag& MeshCell::interior_boundary() const
    {
        return interior_boundary_;
    }

    void MeshCell::mark_to_be_erased()
    {
        erase_ = true;
    }

    bool MeshCell::erase() const
    {
        return erase_;
    }

    void MeshCell::unmark_to_be_erased()
    {
        erase_ = false;
    }

    void MeshCell::set_oga_cell_type(OGA_cell_type_t t)
    {
        OGA_cell_type_ = t;
    }

    //bool MeshCell::is_interior() const
    //{
        //return is_interior_;
    //}

    //bool MeshCell::is_wall() const
    //{
        //return is_wall_;
    //}

    //bool MeshCell::is_dirichlet() const
    //{
        //return is_dirichlet_;
    //}
    
    bool MeshCell::near_boundary() const
    {
        if (!wall_boundary_.empty() || !dirichlet_boundary_.empty() || !empty_boundary_.empty() || !farfield_boundary_.empty())
        {
            return true;
        }

        return false;
    }

    /*bool MeshCell::near_wall() const
    {
        if (!wall_boundary_.empty())
        {
            return true;
        }

        return false;
    }

    bool MeshCell::near_dirichlet() const
    {
        if (!dirichlet_boundary_.empty())
        {
            return true;
        }

        return false;
    }*/

    //const std::vector<Tag>& MeshCell::wall_boundary() const
    const MeshCell::bouarray& MeshCell::wall_boundary() const
    {
        return wall_boundary_;
    }

    //const std::vector<Tag>& MeshCell::dirichlet_boundary() const
    const MeshCell::bouarray& MeshCell::dirichlet_boundary() const
    {
        return dirichlet_boundary_;
    }

    const MeshCell::bouarray& MeshCell::empty_boundary() const
    {
        return empty_boundary_;
    }

    const MeshCell::bouarray& MeshCell::farfield_boundary() const
    {
        return farfield_boundary_;
    }
    
    boundary_t MeshCell::boutype() const
    {
        return boutype_;
    }

    void MeshCell::add_boundary(const Tag& t, boundary_t boutype)
    {
        if (boutype == "wall") {
            wall_boundary_.push_back(t);
        }
        else if (boutype == "dirichlet") {
            dirichlet_boundary_.push_back(t);
        }
        else if (boutype == "empty") {
            empty_boundary_.push_back(t);
        }
        else if (boutype == "farfield") {
            farfield_boundary_.push_back(t);
        }
        else {
            assert(false);
        }
    }

    void MeshCell::set_interior_boundary(const Tag& t)
    {
        interior_boundary_ = t;
    }

    /*void MeshCell::set_partition(bool p)
    {
        is_partition_ = p;
    }
    const bool MeshCell::is_partition() const
    {
        return is_partition_;
    }*/

    void MeshCell::set_boutype(boundary_t btype)
    {
        boutype_ = btype;
    }

    //const std::vector<MeshFace>& MeshCell::face() const
    const MeshCell::facearray& MeshCell::face() const
    {
        return face_;
    }

    //std::vector<MeshFace>& MeshCell::face_p()
    MeshCell::facearray& MeshCell::face_p()
    {
        return face_;
    }

    void MeshCell::make_faces(Shape shape)
    {
        assert(face_.empty());
        assert(!point_.empty());

        if (shape == Shape::quad)
        {
            std::vector<MeshPoint> pts0;
            pts0 = {point_[0], point_[1], point_[2], point_[3]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
        }
        else if (shape == Shape::tri)
        {
            std::vector<MeshPoint> pts0;
            pts0 = {point_[0], point_[1], point_[2]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
        }
        else if (shape == Shape::tet)
        {
            std::vector<MeshPoint> pts0, pts1, pts2, pts3;
            pts0 = {point_[1], point_[2], point_[0]};
            pts1 = {point_[0], point_[3], point_[1]};
            pts2 = {point_[0], point_[2], point_[3]};
            pts3 = {point_[2], point_[1], point_[3]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts1));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts2));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts3));
            assert(!face_.back().face().vertex().empty());
        }
        else if (shape == Shape::pri)
        {
            std::vector<MeshPoint> pts0, pts1, pts2, pts3, pts4;
            //pts0 = {point_[1], point_[0], point_[2]};
            pts0 = {point_[0], point_[1], point_[2]};
            pts1 = {point_[5], point_[4], point_[3]};
            pts2 = {point_[1], point_[0], point_[3], point_[4]};
            pts3 = {point_[3], point_[0], point_[2], point_[5]};
            pts4 = {point_[2], point_[1], point_[4], point_[5]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts1));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts2));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts3));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts4));
            assert(!face_.back().face().vertex().empty());
        }
        else if (shape == Shape::hex)
        {
            std::vector<MeshPoint> pts0, pts1, pts2, pts3, pts4, pts5;
            pts0 = {point_[0], point_[1], point_[2], point_[3]};
            pts1 = {point_[6], point_[5], point_[4], point_[7]};
            pts2 = {point_[2], point_[1], point_[5], point_[6]};
            pts3 = {point_[4], point_[0], point_[3], point_[7]};
            pts4 = {point_[1], point_[0], point_[4], point_[5]};
            pts5 = {point_[3], point_[2], point_[6], point_[7]};
            face_.push_back(MeshFace(pts0));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts1));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts2));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts3));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts4));
            assert(!face_.back().face().vertex().empty());
            face_.push_back(MeshFace(pts5));
            assert(!face_.back().face().vertex().empty());
        }
        else
        {
            assert(false);
        }
    }

    const std::pair<Tag, Tag>& MeshCell::donor() const
    {
        return donor_;
    }

    std::vector<Point> MeshCell::geom_point() const
    {
        std::vector<Point> pts;

        for (const MeshPoint& mp: point_)
            pts.push_back(mp.p());

        return pts;
    }

    const Tag& MeshCell::parent_mesh() const
    {
        return parent_mesh_;
    }

    void MeshCell::remove_all_neighbors()
    {
        //nei_.clear();
        pnei_.clear();
    }

    void MeshCell::remove_parent_cells_of_vertices()
    {
        for (MeshPoint& _p: point_)
        {
            _p.remove_parent_cells();
        }
    }

    void MeshCell::set_point_tag(int i, const Tag& t)
    {
        assert(i >= 0);
        assert(i < point_.size());

        point_[i].set_tag(t);
    }

    /*const bool MeshCell::is_resident() const
      {
      return residency_;
      }*/

    void MeshCell::rotate_points(double ang, int axis, const vec3<double>& rot_axis)
    {
        for (MeshFace& mf: face_)
        {
            mf.rotate(ang, axis, rot_axis);
        }
        for (MeshPoint& point: point_)
        {
            point.rotate_point(ang, axis, rot_axis);
        }

        poly_.rotate_points(ang, axis, rot_axis);
        //polytope_->rotate_points(ang, rot_axis);
    }

    void MeshCell::move_points(const vec3<double>& v)
    {
        for (MeshFace& mf: face_)
        {
            mf.move(v);
        }
        for (MeshPoint& point: point_)
        {
            point.move_point(v);
        }

        poly_.move_points(v);
        //polytope_->move_points(v);
    }

    OGA_cell_type_t MeshCell::oga_cell_type() const
    {
        return OGA_cell_type_;
    }

    void MeshCell::deparent_self_from_vertices()
    {
        for (MeshPoint& p: point_)
            p.remove_parent_cell(tag_);
    }

    void MeshCell::deparent_self_from_faces()
    {
        for (MeshFace& mf: face_)
        {
            mf.remove_parent_cell(tag_);
        }
    }

    void MeshCell::deparent_cell_from_faces(const Tag& tag)
    {
        for (MeshFace& mf: face_)
        {
            mf.remove_parent_cell(tag);
        }
    }

    /*void MeshCell::remove_self_from_neighbors(const Tag& t)
      {
      assert(t() == -1);
      nei_.erase(std::remove(nei_.begin(), nei_.end(), t), nei_.end());
      }*/

    void MeshCell::remove_neighbor(const Tag& t)
    {
        assert(t.isvalid());
        //nei_.erase(std::remove(nei_.begin(), nei_.end(), t), nei_.end());
        //pnei_.erase(std::remove(pnei_.begin(), pnei_.end(), t), pnei_.end());
        pnei_.erase(t);
    }

    std::vector<std::vector<Tag>> MeshCell::face_vertex()
    {
        std::vector<std::vector<Tag>> face_vertex;

        if (point_.size() == 2)
        {
            face_vertex.resize(1);

            face_vertex[0].push_back(point_[0].tag());
            face_vertex[0].push_back(point_[1].tag());
        }
        else if (point_.size() == 3)
        {
            face_vertex.resize(3);

            face_vertex[0].push_back(point_[0].tag());
            face_vertex[0].push_back(point_[1].tag());

            face_vertex[1].push_back(point_[1].tag());
            face_vertex[1].push_back(point_[2].tag());

            face_vertex[2].push_back(point_[2].tag());
            face_vertex[2].push_back(point_[0].tag());
        }
        else if (point_.size() == 4)
        {
            face_vertex.resize(4);

            face_vertex[0].push_back(point_[0].tag());
            face_vertex[0].push_back(point_[1].tag());

            face_vertex[1].push_back(point_[1].tag());
            face_vertex[1].push_back(point_[2].tag());

            face_vertex[2].push_back(point_[2].tag());
            face_vertex[2].push_back(point_[3].tag());

            face_vertex[3].push_back(point_[3].tag());
            face_vertex[3].push_back(point_[0].tag());
        }
        else
        {
            std::cout << "invalid point size in face_vertex()" << std::endl;
            exit(0);
        }

        return face_vertex;
    }

    void MeshCell::set_tag(const Tag& t)
    {
        assert(t.isvalid());
        //if (!first_tag_.isvalid())
        //first_tag_ = t;
        tag_ = t;
    }
    void MeshCell::set_parent_mesh(const Tag& celltag)
    {
        assert(celltag.isvalid());
        //if (!root_parent_mesh_.isvalid())
        //root_parent_mesh_ = celltag;
        parent_mesh_ = celltag;
    }
    const Polyhedron& MeshCell::poly() const
    {
        return poly_;
    }
    /*const Polytope& MeshCell::polytope() const
      {
      return *polytope_;
      }*/
    /*std::shared_ptr<Polytope> MeshCell::polytope() const
      {
      return polytope_;
      }*/
    //const std::vector<MeshPoint>& MeshCell::point() const
    //const std::array<MeshPoint, NPOINT>& MeshCell::point() const
    //{
    //return point_;
    //}
    const MeshCell::pointarray& MeshCell::point() const
    {
        return point_;
    }
    const MeshPoint& MeshCell::point(int i) const
    {
        return point_[i];
    }
    /*const std::vector<Tag>& MeshCell::ipoint() const
      {
      return ipoint_;
      }*/
    /*const std::vector<Tag>& MeshCell::nei() const
      {
      return nei_;
      }*/
    const MeshCell::neiarray& MeshCell::pnei() const
    {
        return pnei_;
    }
    /*const Tag& MeshCell::ipoint(int i) const
      {
      return ipoint_[i];
      }*/
    /*const Tag& MeshCell::nei(int i) const
      {
      return nei_[i];
      }*/
    const Tag& MeshCell::pnei(int i) const
    {
        return pnei_[i];
    }
    const Tag& MeshCell::tag() const
    {
        return tag_;
    }

    MeshCell::MeshCell(const Tag& tag, const Tag& parent_mesh, const std::vector<MeshPoint>& point, boundary_t btype, Shape shape): tag_(tag), parent_mesh_(parent_mesh), OGA_cell_type_(OGA_cell_type_t::undefined), boutype_(btype)
    {
        assert(tag.isvalid());

        //if (btype == boundary_t::wall) {
            //is_wall_ = true;
        //}
        //else if (btype == boundary_t::dirichlet) {
            //is_dirichlet_ = true;
        //}
        //else {
            //is_interior_ = true;
        //}

        point_ = point;
        //assert(point.size() <= NPOINT);
        //std::copy(point.begin(), point.end(), point_.begin());
        //this->point_ = point;
        make_faces(shape);

        std::vector<Point> pts;
        for (const MeshPoint& mp: point)
        {
            pts.push_back(mp.p());
        }

        if (shape == Shape::tet)
        {
            assert(point.size() == 4);
            assert(pts.size() == 4);
        }

        if (pts.size() > 8)
        {
            std::cout << "pts size: " << pts.size() << std::endl;
        }
        assert(pts.size() <= 8);
        poly_ = Polyhedron(pts, shape);
        //poly_.set_vertex(pts);
        //poly_.set_edge();

        //first_tag_ = tag_;
        //root_parent_mesh_ = parent_mesh_;
    }

    void MeshCell::add_vertex(const MeshPoint& p, const Tag& pointtag)
    {
        point_.add(p);
        //assert(npoint_ < NPOINT);
        //point_[npoint_] = p;
        //++npoint_;
        //point_.push_back(p);
    }

    void MeshCell::add_pnei(const Tag& celltag)
    {
        pnei_.push_back(celltag);
    }

    void MeshCell::set_parent_cell_of_vertices()
    {
        for (MeshPoint& _p: point_)
        {
            assert(tag_.isvalid());
            _p.add_parent_cell(tag_);
        }
    }
}
