#include "geom.h"

namespace Common
{
    Polygon create_with_sweep(const Segment& segment, const vec3<double>& v)
    {
        assert(segment.vertex().size() == 2);
        std::vector<Point> vtx;
        std::cout << "pushing vertices" << std::endl;
        vtx.reserve(4);
        vtx.push_back(segment.vertex(0));
        vtx.push_back(segment.vertex(1));
        vtx.push_back(Point(vtx[1].r() + v));
        vtx.push_back(Point(vtx[0].r() + v));
        std::cout << "pushed vertices" << std::endl;
        assert(vtx.size() == 4);
        std::cout << "vvvv: " << vtx[0].r(0) << std::endl;
        std::cout << "vvvv: " << vtx[0].r(1) << std::endl;
        std::cout << "vvvv: " << vtx[1].r(0) << std::endl;
        std::cout << "vvvv: " << vtx[1].r(1) << std::endl;
        std::cout << "vvvv: " << vtx[2].r(0) << std::endl;
        std::cout << "vvvv: " << vtx[2].r(1) << std::endl;
        std::cout << "vvvv: " << vtx[3].r(0) << std::endl;
        std::cout << "vvvv: " << vtx[3].r(1) << std::endl;

        return Polygon(vtx);
    }

    /*bool Polygon::do_intersect(const vec3<double>& v) const
    {
    // any need for this function?
        for (int i=0; i<3; ++i)
        {
            if (r(i) < min_(i) || r(i) > max_(i)) {
                return false;
            }
        }

        Point p_(minx - 5., miny - 5., minz - 5.);
        Point p(r);
        Segment s(p_, p); 

        int counter = 0;
        for (const Edge& f: face_)
        {
            if (f.do_intersect(s)) {
                ++counter;
            }
        }

        if (counter % 2 == 0) {
            return false;
        }

        return true;
    }*/

    const Segment& Polygon::edge(int i) const
    {
        assert(i >= 0);
        assert(i < 2);
        return edge_[i];
    }

    void Polygon::rotate_points(double angle, double axis, const vec3<double>& rot_point)
    {
        RotationMatrix rm;

        for (Point& p: vertex_)
        {
            auto z = p.r() - rot_point;
            auto newz = rm.rotate(angle, axis, z);
            p.set_r(newz + rot_point);
        }

        set_bbox();

        for (Segment& e: edge_)
        {
            e.rotate_points(angle, axis, rot_point);
        }
    }

    void Polygon::move_points(const vec3<double>& v)
    {
        for (Point& _p: vertex_)
        {
            _p.set_r(_p.r(0) + v(0), _p.r(1) + v(1), _p.r(2) + v(2));
        }

        for (Segment& e: edge_)
        {
            e.move_points(v);
        }
    }

    const Polygon::Vertices& Polygon::vertex() const
    {
        return vertex_;
    }

    const Polygon::Edges& Polygon::edge() const
    {
        return edge_;
    }

    void Polygon::set_bbox()
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

    bool Polygon::degenerate() const
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

        int ncop = 0;
        for (int i=0; i<NDIM; ++i)
        {
            if (coplanar(i))
            {
                ++ncop;
                if (ncop == 2)
                {
                    return true;
                }
            }
        }

        return false;
    }

    Polygon::Polygon(const std::vector<Point>& vertex)
    {
        set_vertices(vertex); // because copy constructor cannot be triggered with init list: vertex_(vertex).
        if (!degenerate())
        {
            if (vertex_.size() < 3)
            {
                std::cout << "polygon vertex size: " << vertex_.size() << std::endl;
            }
            assert(vertex_.size() >= 3);
            set_edge();
            set_bbox();
            normal_ = normal();
        }
    }

    void Polygon::set_vertices(const std::vector<Point>& vertex)
    {
        vertex_ = vertex;
    }

    double Polygon::normal(int i) const
    {
        if (normal_) {
            return (*normal_)(i);
        }
        else
        {
            normal();
            return (*normal_)(i);
        }
    }

    vec3<double> Polygon::normal() const
    {
        // 3 points or two edges are enough to calculate normal.
        // No need to calculate area first.

        if (normal_) {
            return *normal_;
        }

        if (vertex_.size() < 3)
        {
            std::cout << "vertex size in normal: " << vertex_.size() << std::endl;
        }
        assert(vertex_.size() >= 3);
 
        vec3<double> a = vertex_[0].r();
        vec3<double> b = vertex_[1].r();
        vec3<double> c = vertex_[2].r();

        auto cr = cross((a-b), (c-b)); 

        if (std::abs(cr(0)) <= ZERO)
        {
            cr.set_x(0.);
        }
        if (std::abs(cr(1)) <= ZERO)
        {
            cr.set_y(0.);
        }
        if (std::abs(cr(2)) <= ZERO)
        {
            cr.set_z(0.);
        }

        normal_ = cr / cr.len();

        return *normal_;
    }

    void Polygon::set_edge()
    {
        assert(edge_.empty());

        // vertex_ must be ready at this point.
        
        if (vertex_.size() == 2)
        {
            edge_.push_back(Segment(vertex_[0], vertex_[1]));
        }
        else if (vertex_.size() == 3)
        {
            edge_.push_back(Segment(vertex_[0], vertex_[1]));
            edge_.push_back(Segment(vertex_[1], vertex_[2]));
            edge_.push_back(Segment(vertex_[2], vertex_[0]));
        }
        else if (vertex_.size() == 4)
        {
            edge_.push_back(Segment(vertex_[0], vertex_[1]));
            edge_.push_back(Segment(vertex_[1], vertex_[2]));
            edge_.push_back(Segment(vertex_[2], vertex_[3]));
            edge_.push_back(Segment(vertex_[3], vertex_[0]));
        }

        assert(vertex_.size() == 2 || vertex_.size() == 3 || vertex_.size() == 4);
    }

    std::vector<double> Polygon::centroid_ab(int i0, int i1) const
    {
        // Taken from https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
        
        std::vector<double> cnt(2);

        cnt[0] = 0.;
        cnt[1] = 0.;

        std::vector<double> m(4);

        for (int i=0; i<vertex_.size()-1; ++i)
        {
            m[0] = vertex_[i].r(i0);
            m[1] = vertex_[i].r(i1);
            m[2] = vertex_[i+1].r(i0);
            m[3] = vertex_[i+1].r(i1);

            double det = m[0] * m[3] - m[1] * m[2];

            cnt[0] += (m[0] + m[2]) * det;
            cnt[1] += (m[1] + m[3]) * det;
        }

        m[0] = vertex_.back().r(i0);
        m[1] = vertex_.back().r(i1);
        m[2] = vertex_[0].r(i0);
        m[3] = vertex_[0].r(i1);

        double det = m[0] * m[3] - m[1] * m[2];

        cnt[0] += (m[0] + m[2]) * det;
        cnt[1] += (m[1] + m[3]) * det;

        assert(!std::isnan(cnt[0]));
        assert(!std::isnan(cnt[1]));

        double oldbef0 = cnt[0];
        double oldbef1 = cnt[1];

        if (vertex_.size() < 3)
        {
            std::cout << "vertex size in centroid ab: " << vertex_.size() << std::endl;
        }
        assert(vertex_.size() >= 3);
        cnt[0] /= (6.*signed_area());       
        cnt[1] /= (6.*signed_area());       

        if (cnt[0] > max_(0))
        {
            std::cout << "cnt[0]: " << cnt[0] << std::endl;
            std::cout << "min_[0]: " << min_(0) << std::endl;
            std::cout << "max_[0]: " << max_(0) << std::endl;
            std::cout << "area: " << signed_area() << std::endl;

            std::cout << "cntbefore[0]: " << oldbef0 << std::endl;
            std::cout << "cntbefore[1]: " << oldbef1 << std::endl;

            std::cout << "r0: " << vertex_[0].r(0) << std::endl;
            std::cout << "r1: " << vertex_[0].r(1) << std::endl;
            std::cout << "r2: " << vertex_[0].r(2) << std::endl;

            std::cout << "r0: " << vertex_[1].r(0) << std::endl;
            std::cout << "r1: " << vertex_[1].r(1) << std::endl;
            std::cout << "r2: " << vertex_[1].r(2) << std::endl;

            std::cout << "r0: " << vertex_[2].r(0) << std::endl;
            std::cout << "r1: " << vertex_[2].r(1) << std::endl;
            std::cout << "r2: " << vertex_[2].r(2) << std::endl;

            std::cout << "r0: " << vertex_[3].r(0) << std::endl;
            std::cout << "r1: " << vertex_[3].r(1) << std::endl;
            std::cout << "r2: " << vertex_[3].r(2) << std::endl;
        }
        assert(cnt[0] <= max_(0));

        assert(signed_area() != 0.);
        assert(!std::isnan(signed_area()));

        assert(!std::isnan(cnt[0]));
        assert(!std::isnan(cnt[1]));

        return cnt;
    }

    vec3<double> Polygon::centroid() const
    {
        // https://stackoverflow.com/a/2360507/1128551


        // for now I am using simple arithmetic average.
        // should work good for simple shapes such as triangle and quad.

        double cntx = 0.;
        double cnty = 0.;
        double cntz = 0.;
        for (const auto& v: vertex_)
        {
            cntx += v.r(0);
            cnty += v.r(1);
            cntz += v.r(2);
        }

        vec3<double> cnt;
        cnt.set_x(cntx / vertex_.size());
        cnt.set_y(cnty / vertex_.size());
        cnt.set_z(cntz / vertex_.size());

        /*std::vector<double> cnt_xy = Polygon::centroid_ab(0, 1);
        std::vector<double> cnt_xz = Polygon::centroid_ab(0, 2);

        vec3<double> cnt;
        cnt.set_x(cnt_xy[0]);
        cnt.set_y(cnt_xy[1]);
        cnt.set_z(cnt_xz[1]);*/

        assert(!std::isnan(cnt(0)));
        assert(!std::isnan(cnt(1)));
        assert(!std::isnan(cnt(2)));

        //if (cnt(1) > max_(1))
        //{
            //std::cout << "cnt1: " << cnt(1) << std::endl;
            //std::cout << "maxx: " << max_(1) << std::endl;
        //}

        if ((cnt(0) - max_(0)) >= ZERO)
        {
            for (int i=0; i<vertex_.size(); ++i)
            {
                std::cout << "aaa: " << vertex_[i].r(0) << std::endl;
                std::cout << "aaa: " << vertex_[i].r(1) << std::endl;
                std::cout << "aaa: " << vertex_[i].r(2) << std::endl;
            }
            std::cout << "cnt0: " << cnt(0) << std::endl;
            std::cout << "maxx: " << max_(0) << std::endl;
        }
        assert((cnt(0) - max_(0)) < ZERO);
        assert((cnt(1) - max_(1)) < ZERO);
        assert((cnt(2) - max_(2)) < ZERO);

        assert((-cnt(0) + min_(0)) < ZERO);
        assert((-cnt(1) + min_(1)) < ZERO);
        assert((-cnt(2) + min_(2)) < ZERO);

        return cnt;
    }

    double Polygon::signed_area() const
    {
        assert(!vertex_.empty());
        int n = vertex_.size();
        double area = 0.;
        double an, ax, ay, az; // abs value of normal and its coords
        int coord;           // coord to ignore: 1=x, 2=y, 3=z
        int i, j, k;         // loop indices

        if (vertex_.size() < 3)
        {
            std::cout << "vertex size in signed area: " << vertex_.size() << std::endl;
        }
        assert(vertex_.size() >= 3);
        auto N = normal();

        // select largest abs coordinate to ignore for projection
        ax = (N(0)>0 ? N(0) : -N(0));    // abs x-coord
        ay = (N(1)>0 ? N(1) : -N(1));    // abs y-coord
        az = (N(2)>0 ? N(2) : -N(2));    // abs z-coord

        coord = 3;                    // ignore z-coord
        if (ax > ay) {
            if (ax > az) coord = 1;   // ignore x-coord
        }
        else if (ay > az) coord = 2;  // ignore y-coord

        // compute area of the 2D projection
        assert(vertex_.size() >= 3);
        switch (coord) {
            case 1:
                for (i=1, j=2, k=0; i<n; i++, j++, k++)
                {
                    assert(i <= 3);
                    //assert(j <= 3);
                    assert(k <= 3);
                    area += (vertex_[i].r(1) * (vertex_[j==n ? 0 : j].r(2) - vertex_[k].r(2)));
                }
                break;
            case 2:
                for (i=1, j=2, k=0; i<n; i++, j++, k++)
                {
                    assert(i <= 3);
                    //assert(j <= 3);
                    assert(k <= 3);
                    area += (vertex_[i].r(2) * (vertex_[j==n ? 0 : j].r(0) - vertex_[k].r(0)));
                }
                break;
            case 3:
                for (i=1, j=2, k=0; i<n; i++, j++, k++)
                {
                    assert(i <= 3);
                    //assert(j <= 3);
                    assert(k <= 3);
                    area += (vertex_[i].r(0) * (vertex_[j==n ? 0 : j].r(1) - vertex_[k].r(1)));
                }
                break;
        }
        switch (coord) {    // wrap-around term
            case 1:
                area += (vertex_[0].r(1) * (vertex_[1].r(2) - vertex_.back().r(2)));
                break;
            case 2:
                area += (vertex_[0].r(2) * (vertex_[1].r(0) - vertex_.back().r(0)));
                break;
            case 3:
                area += (vertex_[0].r(0) * (vertex_[1].r(1) - vertex_.back().r(1)));
                break;
        }

        // scale to get area before projection
        an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
        switch (coord) {
            case 1:
                area *= (an / (2 * N(0)));
                break;
            case 2:
                area *= (an / (2 * N(1)));
                break;
            case 3:
                area *= (an / (2 * N(2)));
        }
        signed_area_ = -area;
        return *signed_area_;
    }

    /*double Polygon::signed_area() const
    {
        // Adapted from http://mathworld.wolfram.com/PolygonArea.html.
        //
        if (signed_area_) {
            return *signed_area_;
        }

        assert(vertex_.size() > 0);
        assert(vertex_.size() < 5);

        double sa = 0.;
        //signed_area_(0.);
        //arma::Mat<double>::fixed<2, 2> m;
        std::vector<double> m(4);

        for (int i=0; i<vertex_.size()-1; ++i)
        {
            m[0] = vertex_[i].r(0);
            m[1] = vertex_[i].r(1);
            m[2] = vertex_[i+1].r(0);
            m[3] = vertex_[i+1].r(1);

            double det = m[0] * m[3] - m[1] * m[2];

            //signed_area += arma::det(m);
            sa += det;
        }

        m[0] = vertex_.back().r(0);
        m[1] = vertex_.back().r(1);
        m[2] = vertex_[0].r(0);
        m[3] = vertex_[0].r(1);

        //signed_area += arma::det(m);
        double det = m[0] * m[3] - m[1] * m[2];
        sa += det;
        sa *= 0.5;

        signed_area_ = sa;

        return *signed_area_;
    }*/

    double magcross2D(const vec3<double>& a, const vec3<double>& b)
    {
        return std::abs(a(0) * b(1) - a(1) * b(0));
    }

    //vec3<double> normalize(const vec3<double>& a)
    //{
        //double mag = std::sqrt(a(0)*a(0) + a(1)*a(1));
        //vec3<double> b;
        //b.set(a(0)/mag, b(0)/mag);
        //return b;
    //}

    const Point& Polygon::vertex(int i) const
    {
        return vertex_[i];
    }

    bool Polygon::do_intersect(const Segment& s) const
    {
        // https://math.stackexchange.com/questions/47594/plane-intersecting-line-segment

        vec3<double> p0 = s.vertex(0).r();
        vec3<double> p1 = s.vertex(1).r();
        vec3<double> v0 = vertex_[0].r();
        if (vertex_.size() < 3)
        {
            std::cout << "vertex size in signed area: " << vertex_.size() << std::endl;
        }
        auto n = normal();

        auto u = p0 - v0;
        auto v = p1 - v0;

        double uu = dotp(n, u);
        double vv = dotp(n, v);

        double res = uu * vv;

        //std::cout << "p0: " << p0(0) << std::endl;
        //std::cout << "p0: " << p0(1) << std::endl;
        //std::cout << "p0: " << p0(2) << std::endl;

        //std::cout << "p1: " << p1(0) << std::endl;
        //std::cout << "p1: " << p1(1) << std::endl;
        //std::cout << "p1: " << p1(2) << std::endl;

        //std::cout << "n: " << n(0) << std::endl;
        //std::cout << "n: " << n(1) << std::endl;
        //std::cout << "n: " << n(2) << std::endl;

        //std::cout << "u: " << u(0) << std::endl;
        //std::cout << "u: " << u(1) << std::endl;
        //std::cout << "u: " << u(2) << std::endl;

        //std::cout << "v: " << v(0) << std::endl;
        //std::cout << "v: " << v(1) << std::endl;
        //std::cout << "v: " << v(2) << std::endl;

        //std::cout << "uu: " << uu << std::endl;
        //std::cout << "vv: " << vv << std::endl;

        //std::cout << "res: " << res << std::endl;

        if (res <= ZERO) {
            return true;
        }

        return false;
    }

    /*bool Polygon::do_intersect(const Segment& s) const
    {
        // based on http://geomalgorithms.com/a05-_intersect-1.html

        vec3<double> p0 = s.vertex(0).r();
        vec3<double> p1 = s.vertex(1).r();
        vec3<double> v0 = vertex_[0].r();

        auto w = v0 - p0;
        auto u = p1 - p0;
        auto n = normal();

        // first check if the polygon and the segment are parallel.
        if (std::abs(dotp(n, u)) <= ZERO) {

            std::cout << "segment and polygron are parallel" << std::endl;
            return false;
        }

        std::cout << "segment and polygron are not parallel" << std::endl;

        double si = -dotp(n, w) / dotp(n, u);

            std::cout << "n: " << n(0) << std::endl;
            std::cout << "n: " << n(1) << std::endl;
            std::cout << "n: " << n(2) << std::endl;

            std::cout << "u: " << u(0) << std::endl;
            std::cout << "u: " << u(1) << std::endl;
            std::cout << "u: " << u(2) << std::endl;

            std::cout << "p0: " << p0(0) << std::endl;
            std::cout << "p0: " << p0(1) << std::endl;
            std::cout << "p0: " << p0(2) << std::endl;

            std::cout << "p1: " << p1(0) << std::endl;
            std::cout << "p1: " << p1(1) << std::endl;
            std::cout << "p1: " << p1(2) << std::endl;

            std::cout << "v0: " << v0(0) << std::endl;
            std::cout << "v0: " << v0(1) << std::endl;
            std::cout << "v0: " << v0(2) << std::endl;

            std::cout << "w: " << w(0) << std::endl;
            std::cout << "w: " << w(1) << std::endl;
            std::cout << "w: " << w(2) << std::endl;

            std::cout << "si " << si << std::endl;

        if (si < 0 || si > 1) {
            return false;
        }

        return true;
    }*/

    //std::vector<CGAL_Point> Polygon::cgalpoint() const
    //{
        //std::vector<CGAL_Point> pts;
        //for (const Point& p: vertex_)
        //{
            //pts.push_back(CGAL_Point(p.r(0), p.r(1)));
        //}
//
        //return pts;
    //}

    //Polygon_2 Polygon::cgalpolygon() const
    //{
        //auto pts = this->cgalpoint();
        //return Polygon_2(pts.begin(), pts.end());
    //}

}
