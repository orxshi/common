#include "geom.h"

namespace Common
{
    Polyhedron create_with_sweep(const Polygon& polygon, const vec3<double>& v)
    {
        if (polygon.vertex().size() == 3)
        {
            std::vector<Point> vtx;
            vtx.reserve(6);
            vtx.push_back(polygon.vertex(0));
            vtx.push_back(polygon.vertex(1));
            vtx.push_back(polygon.vertex(2));
            vtx.push_back(Point(vtx[0].r() + v));
            vtx.push_back(Point(vtx[1].r() + v));
            vtx.push_back(Point(vtx[2].r() + v));

            return Polyhedron(vtx, Shape::pri);
        }
        else if (polygon.vertex().size() == 4)
        {
            std::vector<Point> vtx;
            vtx.reserve(8);
            vtx.push_back(polygon.vertex(0));
            vtx.push_back(polygon.vertex(1));
            vtx.push_back(polygon.vertex(2));
            vtx.push_back(polygon.vertex(3));
            vtx.push_back(Point(vtx[0].r() + v));
            vtx.push_back(Point(vtx[1].r() + v));
            vtx.push_back(Point(vtx[2].r() + v));
            vtx.push_back(Point(vtx[3].r() + v));

            return Polyhedron(vtx, Shape::hex);
        }
        else
        {
            assert(false);
        }
    }

    size_t Polyhedron::mem() const
    {
        size_t size = 0;
        size += vertex_.mem();
        size += face_.mem();
        size += sizeof(double);
        size += sizeof(int);
        size += sizeof(double) * 6;

        return size;
    }

    Polyhedron::Polyhedron(const std::vector<Point>& vertex, Shape shape): shape_(shape)
    {
        if (vertex.size() > 8)
        {
            std::cout << "vertex size: " << vertex.size() << std::endl;
            std::cout << "shape: " << static_cast<int>(shape) << std::endl;
        }
        assert(vertex.size() <= 8);
        if (shape_ == Shape::tet)
        {
            assert(vertex.size() == 4);
        }
        set_vertices(vertex); // because copy constructor cannot be triggered with init list: vertex_(vertex).
        set_faces();
        set_bbox();
        assert(!face_.empty());
    }

    Shape Polyhedron::shape() const
    {
        return shape_;
    }

    void Polyhedron::set_min(vec3<double> min)
    {
        min_ = min;
    }

    void Polyhedron::set_max(vec3<double> max)
    {
        max_ = max;
    }

    void Polyhedron::set_vertices_from_bbox()
    {
        assert(vertex_.size() == 0);

        vertex_.reserve(8);

        vertex_[0].set_r(max(0), max(1), min(2));
        vertex_[1].set_r(min(0), max(1), min(2));
        vertex_[2].set_r(min(0), min(1), min(2));
        vertex_[3].set_r(max(0), min(1), min(2));

        vertex_[4].set_r(max(0), max(1), max(2));
        vertex_[5].set_r(min(0), max(1), max(2));
        vertex_[6].set_r(min(0), min(1), max(2));
        vertex_[7].set_r(max(0), min(1), max(2));
    }

    void Polyhedron::set_bbox(vec3<double> min, vec3<double> max)
    {
        min_ = min;
        max_ = max;
    }

    void Polyhedron::set_bbox(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
    {
        min_.set(xmin, ymin, zmin);
        max_.set(xmax, ymax, zmax);
    }

    void Polyhedron::set_shape(Shape shape)
    {
        shape_ = shape;
    }

    void Polyhedron::set_bbox()
    {
        min_.set_x(BIG_POS_NUM);
        min_.set_y(BIG_POS_NUM);
        min_.set_z(BIG_POS_NUM);
        max_.set_x(BIG_NEG_NUM);
        max_.set_y(BIG_NEG_NUM);
        max_.set_z(BIG_NEG_NUM);

        for (int i=0; i<vertex_.size(); ++i)
        {
            min_.set_x(std::min(min_(0), vertex_[i].r(0)));
            min_.set_y(std::min(min_(1), vertex_[i].r(1)));
            min_.set_z(std::min(min_(2), vertex_[i].r(2)));

            max_.set_x(std::max(max_(0), vertex_[i].r(0)));
            max_.set_y(std::max(max_(1), vertex_[i].r(1)));
            max_.set_z(std::max(max_(2), vertex_[i].r(2)));
        }
    }

    void Polyhedron::set_vertices(const std::vector<Point>& vertex)
    {
        if (vertex.size() > 8)
        {
            std::cout << "vertex size: " << vertex.size() << std::endl;
        }
        assert(vertex.size() <= 8);
        vertex_ = vertex;
        if (shape_ == Shape::tet)
        {
            assert(vertex.size() == 4);
            assert(vertex_.size() == 4);
        }
    }

    const Polygon& Polyhedron::face(int i) const
    {
        assert(i >= 0);
        assert(i < 6);
        return face_[i];
    }

    void Polyhedron::rotate_points(double angle, double axis, const vec3<double>& rot_point)
    {
        RotationMatrix rm;

        for (Point& p: vertex_)
        {
            auto z = p.r() - rot_point;
            auto newz = rm.rotate(angle, axis, z);
            p.set_r(newz + rot_point);
        }

        set_bbox();
        for (Polygon& f: face_)
        {
            f.rotate_points(angle, axis, rot_point);
        }
    }

    const vec3<double>& Polyhedron::min() const
    {
        //std::cout << "returning hedron min" << std::endl;
        //std::cout << "min0: " << min_(0) << std::endl;
        //std::cout << "min1: " << min_(1) << std::endl;
        //std::cout << "min2: " << min_(2) << std::endl;
        return min_;
    }
    
    const vec3<double>& Polyhedron::max() const
    {
        return max_;
    }

    double Polyhedron::min(int i) const
    {
        assert(i < 3);
        assert(i >= 0);
        //std::cout << "returning hedron min i" << std::endl;
        //std::cout << "minn0: " << min_(0) << std::endl;
        //std::cout << "minn1: " << min_(1) << std::endl;
        //std::cout << "minn2: " << min_(2) << std::endl;
        return min_(i);
    }

    double Polyhedron::max(int i) const
    {
        assert(i < 3);
        assert(i >= 0);
        return max_(i);
    }

    void Polyhedron::move_points(const vec3<double>& v)
    {
        for (Point& _p: vertex_)
        {
            _p.set_r(_p.r(0) + v(0), _p.r(1) + v(1), _p.r(2) + v(2));
        }

        set_bbox();

        for (Polygon& f: face_)
        {
            f.move_points(v);
        }
    }

    const Polyhedron::Vertices& Polyhedron::vertices() const
    {
        return vertex_;
    }

    const Polyhedron::Faces& Polyhedron::faces() const
    {
        return face_;
    }

    //double Polyhedron::volume_triface(std::vector<Polygon>::const_iterator begin, std::vector<Polygon>::const_iterator end) const
    double Polyhedron::volume_triface(const Polygon* begin, const Polygon* end) const
    {
        // no need to call directly.
        // use Polyhedron::volume() directly.

        double sum = 0.;

        //std::cout << "****************" << std::endl;
        //std::cout << "nface: " << end - begin << std::endl;
        //std::cout << "----------------" << std::endl;

        for (auto f = begin; f != end; ++f)
        {
            if (f->vertex().size() < 3)
            {
                std::cout << "vertex size in hedron volume triface: " << f->vertex().size() << std::endl;
            }
            assert(f->vertex().size() >= 3);
            assert(!std::isnan(f->signed_area()));
            assert(f->signed_area() != 0.);

            assert(!std::isnan(f->centroid()(0)));
            assert(!std::isnan(f->centroid()(1)));
            assert(!std::isnan(f->centroid()(2)));

            assert(!std::isnan(f->normal()(0)));
            assert(!std::isnan(f->normal()(1)));
            assert(!std::isnan(f->normal()(2)));

            assert(f->signed_area() != 0.);

            //assert(f->centroid()(0) != 0.);
            //assert(f->centroid()(1) != 0.);
            //assert(f->centroid()(2) != 0.);

            //assert(f->normal()(0) != 0.);
            //assert(f->normal()(1) != 0.);
            //assert(f->normal()(2) != 0.);

            //if (dotp(f->centroid(), f->normal()) == 0.)
            //{
                //std::cout << "cnt 0: " << f->centroid()(0) << std::endl;
                //std::cout << "cnt 1: " << f->centroid()(1) << std::endl;
                //std::cout << "cnt 2: " << f->centroid()(2) << std::endl;

                //std::cout << "norm 0: " << f->normal()(0) << std::endl;
                //std::cout << "norm 1: " << f->normal()(1) << std::endl;
                //std::cout << "norm 2: " << f->normal()(2) << std::endl;

                //std::cout << "area: " << std::abs(f->signed_area()) << std::endl;
            //}
            //assert(dotp(f->centroid(), f->normal()) != 0.);
            
            sum += dotp(f->centroid(), f->normal()) * std::abs(f->signed_area());

            //std::cout << "tsum: " << dotp(f->centroid(), f->normal()) * std::abs(f->signed_area()) << std::endl;
            //std::cout << "----------------" << std::endl;
        }

        assert(!std::isnan(sum));
        //assert(sum != 0.);

        double volume = sum / 3.;

        return volume;
    }

    vec3<double> Polyhedron::centroid() const
    {
        // this assumes that only vertices have weight.

        double cnt_x = 0.;
        double cnt_y = 0.;
        double cnt_z = 0.;

        for (const auto& p: vertex_)
        {
            cnt_x += p.r(0);
            cnt_y += p.r(1);
            cnt_z += p.r(2);
        }

        cnt_x /= vertex_.size();
        cnt_y /= vertex_.size();
        cnt_z /= vertex_.size();

        vec3<double> cnt;
        cnt.set_x(cnt_x);
        cnt.set_y(cnt_y);
        cnt.set_z(cnt_z);

        assert(cnt(0) >= min_(0));
        assert(cnt(1) >= min_(1));
        assert(cnt(2) >= min_(2));

        assert(cnt(0) <= max_(0));
        assert(cnt(1) <= max_(1));
        assert(cnt(2) <= max_(2));

        return cnt;
    }

    /*vec3<double> Polyhedron::centroid() const
    {
        // based on http://wwwf.imperial.ac.uk/~rn/centroid.pdf

        double cnt_x = 0.;
        double cnt_y = 0.;
        double cnt_z = 0.;

        for (auto f = face_.begin(); f != face_.end(); ++f)
        //for (auto f = begin; f != end; ++f)
        {
            assert(!std::isnan(f->vertex(0).r(0)));
            assert(!std::isnan(f->vertex(0).r(1)));
            assert(!std::isnan(f->vertex(0).r(2)));

            assert(!std::isnan(f->vertex(1).r(0)));
            assert(!std::isnan(f->vertex(1).r(1)));
            assert(!std::isnan(f->vertex(1).r(2)));

            assert(!std::isnan(f->vertex(2).r(0)));
            assert(!std::isnan(f->vertex(2).r(1)));
            assert(!std::isnan(f->vertex(2).r(2)));

            double cx = (f->vertex(0).r(0) + f->vertex(1).r(0) + f->vertex(2).r(0)) / 3.;
            double cy = (f->vertex(0).r(1) + f->vertex(1).r(1) + f->vertex(2).r(1)) / 3.;
            double cz = (f->vertex(0).r(2) + f->vertex(1).r(2) + f->vertex(2).r(2)) / 3.;

            assert(!std::isnan(cx));
            assert(!std::isnan(cy));
            assert(!std::isnan(cz));

            vec3<double> n = f->normal();
            assert(!std::isnan(n(0)));
            assert(!std::isnan(n(1)));
            assert(!std::isnan(n(2)));

            cnt_x += (cx * cx) * n(0);
            cnt_y += (cy * cy) * n(1);
            cnt_z += (cz * cz) * n(2);
        }

        vec3<double> cnt;
        cnt.set_x(cnt_x);
        cnt.set_y(cnt_y);
        cnt.set_z(cnt_z);

        assert(!std::isnan(cnt_x));
        assert(!std::isnan(cnt_y));
        assert(!std::isnan(cnt_z));

        assert(volume() != 0.);
        assert(!std::isnan(volume()));
        cnt = cnt / (2. * volume());

        if (cnt(0) < min_(0))
        {
            std::cout << "cnt(0): " << cnt(0) << std::endl;
            std::cout << "cnt(1): " << cnt(1) << std::endl;
            std::cout << "cnt(2): " << cnt(2) << std::endl;

            std::cout << "min(0): " << min_(0) << std::endl;
            std::cout << "min(1): " << min_(1) << std::endl;
            std::cout << "min(2): " << min_(2) << std::endl;

            std::cout << "max(0): " << max_(0) << std::endl;
            std::cout << "max(1): " << max_(1) << std::endl;
            std::cout << "max(2): " << max_(2) << std::endl;

            std::cout << "vol: " << volume() << std::endl;

            std::cout << "shape: " << static_cast<int>(shape_) << std::endl;

            for (const Point& p: vertex_)
            {
                std::cout << "r: " << p.r(0) << std::endl;
                std::cout << "r: " << p.r(1) << std::endl;
                std::cout << "r: " << p.r(2) << std::endl;
                std::cout << std::endl;
            }

            for (const Polygon& f: face_)
            {
                for (const Point& p: f.vertex())
                {
                    std::cout << "vertex: " << p.r(0) << " " << p.r(1) << " " << p.r(2) << std::endl;
                }
                std::cout << "norm: " << f.normal()(0) << " " << f.normal()(1) << " " << f.normal()(2) << std::endl;
                std::cout << "signedarea: " << f.signed_area() << std::endl;
            }
        }

        assert(cnt(0) >= min_(0));
        assert(cnt(1) >= min_(1));
        assert(cnt(2) >= min_(2));

        assert(cnt(0) <= max_(0));
        assert(cnt(1) <= max_(1));
        assert(cnt(2) <= max_(2));

        return cnt;
    }*/

    bool Polyhedron::degenerate() const
    {
        auto coplanar = [&](int j)
        {
            bool cop = true;
            for (int i=0; i<vertex_.size()-1; ++i)
            {
                if (std::abs(vertex_[i].r(j) - vertex_[i+1].r(j)) > ZERO)
                {
                    cop = false;
                    break;
                }
            }
            return cop;
        };

        for (int i=0; i<NDIM; ++i)
        {
            if (coplanar(i))
            {
                return true;
            }
        }

        return false;
    }

    double Polyhedron::volume() const
    {
        // based on https://cfd.ku.edu/papers/1999-AIAAJ.pdf

        if (volume_) {
            return *volume_;
        }

        double volume;

        if (shape_ == Shape::tet)
        {
            if (degenerate())
            {
                return 0.;
            }
            for (const auto& ff: face_)
            {
                if (ff.vertex().size() < 3)
                {
                    std::cout << "vertex size in hedron volume triface: " << ff.vertex().size() << std::endl;
                }
                assert(ff.vertex().size() >= 3);
            }

            //volume = volume_triface(face_.begin(), face_.end());
            volume = volume_triface(face_.pbegin(), face_.pend());
        }
        else if (shape_ == Shape::hex)
        {
            if (degenerate())
            {
                return 0.;
            }
            // break each face of hex into 2 triangles.
            // add triangles to container.

            std::array<Polygon, 12> faces;

            assert(!vertex_.empty());
            assert(!face_.empty());

            int i = 0;
            for (const Polygon& f: face_)
            {
                faces[i]   = Polygon(std::vector<Point>{f.vertex(0), f.vertex(1), f.vertex(2)});
                assert(faces[i].vertex().size() >= 3);
                faces[i+1] = Polygon(std::vector<Point>{f.vertex(0), f.vertex(2), f.vertex(3)});
                assert(faces[i+1].vertex().size() >= 3);

                i = i + 2;
            }

            for (const auto& ff: faces)
            {
                if (ff.vertex().size() < 3)
                {
                    std::cout << "vertex size in hedron volume triface: " << ff.vertex().size() << std::endl;
                }
                assert(ff.vertex().size() >= 3);
            }

            volume = volume_triface(faces.begin(), faces.end());
        }
        else if (shape_ == Shape::pri)
        {
            if (degenerate())
            {
                return 0.;
            }
            // break each quad face of prism into 2 triangles: 2 tri + 2 * 3 tri = 8 tri
            // add triangles to container.

            std::array<Polygon, 8> faces;

            int i = 0;
            for (const Polygon& f: face_)
            {
                if (f.vertex().size() == 3)
                {
                    faces[i] = f;
                    ++i;
                }
                else
                {
                    faces[i]   = Polygon(std::vector<Point>{f.vertex(0), f.vertex(1), f.vertex(2)});
                    faces[i+1] = Polygon(std::vector<Point>{f.vertex(0), f.vertex(2), f.vertex(3)});
                    i = i + 2;
                }
            }

            for (const auto& ff: faces)
            {
                if (ff.vertex().size() < 3)
                {
                    std::cout << "vertex size in hedron volume triface: " << ff.vertex().size() << std::endl;
                }
                assert(ff.vertex().size() >= 3);
            }

            volume = volume_triface(faces.begin(), faces.end());
        }

        return std::abs(volume);
    }

    void Polyhedron::set_faces()
    {
        if (vertex_.empty()) {
            return;
        }

        if (shape_ == Shape::quad)
        {
            if (vertex_.size() != 4)
            {
                std::cout << "vertex size: " << vertex_.size() << std::endl;
            }
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[1], vertex_[2], vertex_[3]}));
        }
        else if (shape_ == Shape::tri)
        {
            if (vertex_.size() != 3)
            {
                std::cout << "vertex size: " << vertex_.size() << std::endl;
            }
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[1], vertex_[2]}));
        }
        else if (shape_ == Shape::tet)
        {
            if (vertex_.size() != 4)
            {
                std::cout << "vertex size: " << vertex_.size() << std::endl;
            }
            assert(vertex_.size() == 4);
            face_.push_back(Polygon(std::vector<Point>{vertex_[1], vertex_[2], vertex_[0]})); // back face
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[3], vertex_[1]})); // front left
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[2], vertex_[3]})); // front right
            face_.push_back(Polygon(std::vector<Point>{vertex_[2], vertex_[1], vertex_[3]})); // bottom
        }
        else if (shape_ == Shape::pri)
        {
            assert(vertex_.size() == 6);
            //face_.push_back(Polygon(std::vector<Point>{vertex_[1], vertex_[0], vertex_[2]})); // back face
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[1], vertex_[2]})); // back face
            face_.push_back(Polygon(std::vector<Point>{vertex_[5], vertex_[4], vertex_[3]})); // front face
            face_.push_back(Polygon(std::vector<Point>{vertex_[1], vertex_[0], vertex_[3], vertex_[4]})); // left face
            face_.push_back(Polygon(std::vector<Point>{vertex_[3], vertex_[0], vertex_[2], vertex_[5]})); // oblique face
            face_.push_back(Polygon(std::vector<Point>{vertex_[2], vertex_[1], vertex_[4], vertex_[5]})); // bottom face
        }
        else if (shape_ == Shape::hex)
        {
            assert(vertex_.size() == 8);
            face_.push_back(Polygon(std::vector<Point>{vertex_[0], vertex_[1], vertex_[2], vertex_[3]})); // back face
            face_.push_back(Polygon(std::vector<Point>{vertex_[6], vertex_[5], vertex_[4], vertex_[7]})); // front face
            face_.push_back(Polygon(std::vector<Point>{vertex_[2], vertex_[1], vertex_[5], vertex_[6]})); // left face
            face_.push_back(Polygon(std::vector<Point>{vertex_[4], vertex_[0], vertex_[3], vertex_[7]})); // right face
            face_.push_back(Polygon(std::vector<Point>{vertex_[1], vertex_[0], vertex_[4], vertex_[5]})); // top face
            face_.push_back(Polygon(std::vector<Point>{vertex_[3], vertex_[2], vertex_[6], vertex_[7]})); // bottom face
        }
        else 
        {
            assert(false);
        }
        assert(!face_.empty());
    }

    bool Polyhedron::do_intersect(const vec3<double>& r) const
    {
        for (int i=0; i<3; ++i)
        {
            if (min_(i) - r(i) > ZERO) {
                return false;
            }

            if (r(i) - max_(i) > ZERO) {
                return false;
            }

            //if (r(i) < min_(i) || r(i) > max_(i)) {
                //return false;
            //}
        }

        //std::cout << "r is inside aabb" << std::endl;

        Point p_(min_(0) - 5., min_(1) - 5., min_(2) - 5.);
        Point p(r);
        Segment s(p_, p); 
        assert(s.vertex(0).r(0) == (min_(0) - 5.));
        assert(s.vertex(0).r(1) == (min_(1) - 5.));
        assert(s.vertex(0).r(2) == (min_(2) - 5.));

        int counter = 0;
        assert(!face_.empty());
        for (const Polygon& f: face_)
        {
            //std::cout << "checking inter" << std::endl;
            //std::cout << "face vertex: " << f.vertex(0).r(0) << " " << f.vertex(0).r(1) << " " << f.vertex(0).r(2) << std::endl;
            //std::cout << "face vertex: " << f.vertex(1).r(0) << " " << f.vertex(1).r(1) << " " << f.vertex(1).r(2) << std::endl;
            //std::cout << "face vertex: " << f.vertex(2).r(0) << " " << f.vertex(2).r(1) << " " << f.vertex(2).r(2) << std::endl;
            //std::cout << "face vertex: " << f.vertex(3).r(0) << " " << f.vertex(3).r(1) << " " << f.vertex(3).r(2) << std::endl;

            //std::cout << "seg vertex: " << s.vertex(0).r(0) << " " << s.vertex(0).r(1) << " " << s.vertex(0).r(2) << std::endl;
            //std::cout << "seg vertex: " << s.vertex(1).r(0) << " " << s.vertex(1).r(1) << " " << s.vertex(1).r(2) << std::endl;

            //std::cout << "min_: " << min_(0) << " " << min_(1) << " " << min_(2) << std::endl;
            //std::cout << "min: " << min(0) << " " << min(1) << " " << min(2) << std::endl;

            //std::cout << "min_ 5: " << min_(0) - 5. << " " << min_(1) - 5. << " " << min_(2) - 5. << std::endl;

            if (f.do_intersect(s)) {
                //std::cout << "yesinter" << std::endl;
                ++counter;
            }
            //else
            //{
                //std::cout << "nointer" << std::endl;
            //}
        }

        //std::cout << "ninter: " << counter << std::endl;

        if (counter % 2 == 0) {
            return false;
        }

        return true;
    }

    bool Polyhedron::do_intersect(const Segment& s) const
    {
        for (const Polygon& f: face_)
        {
            if (f.do_intersect(s)) {
                return true;
            }
        }

        for (int i=0; i<NDIM; ++i)
        {
            if (s.min(i) < min_(i)) {
                return false;
            }

            if (s.max(i) > max_(i)) {
                return false;
            }
        }        

        return true;
    }
}
