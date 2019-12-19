#include "geom.h"

namespace Common
{
    size_t Point::mem() const
    {
        size_t size = 0;

        size += sizeof(double) * 3;

        return size;
    }

    //CGAL_Point Point::cgal_point() const
    //{
        //return CGAL_Point(r_(0), r_(1));
    //}

    bool Point::operator==(const Point& other) const
    {
        return r_ == other.r();
    }

    Point::Point(const Point& p)
    {
        r_ = p.r_;
        assert(r_(0) == p.r_(0));
    }

    Point& Point::operator=(const Point& p)
    {
        r_ = p.r_;
        assert(r_(0) == p.r_(0));
        return *this;
    }

    const vec3<double>& Point::r() const
    {
        return r_;
    }

    double Point::r(int i) const
    {
        assert(i >= 0);
        assert(i < 3);
        return r_(i);
    }

    void Point::set_r(double x, double y, double z)
    {
        r_.set(x, y, z);
    }

    void Point::set_r(const vec3<double>& _r)
    {
        r_.set(_r(0), _r(1), _r(2));
    }

    Point::Point(double x, double y, double z): r_(x, y, z)
    {
    }

    Point::Point(const vec3<double>& r): r_(r)
    {
    }
}
