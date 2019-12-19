#ifndef COMMON_GEOMETRIC_SHAPE_H
#define	COMMON_GEOMETRIC_SHAPE_H

#include <vector>
#include <algorithm>
#include "commoncore.h"
#include "vec3.h"
#include "array.h"
#include "tag.h"
#include "rotationmatrix.h"

namespace Common
{
    class Point
    {
        vec3<double> r_;
        friend class boost::serialization::access;

        public:

        Point(): Point(0., 0., 0.) {};
        Point(double x, double y, double z);
        Point(const vec3<double>& r);
        Point(const Point& p);

        size_t mem() const;
        const vec3<double>& r() const;
        double r(int i) const;
        void set_r(double x, double y, double z);
        void set_r(const vec3<double>& _r);
        bool operator==(const Point& other) const;
        Point& operator=(const Point& p);
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & r_;
        }
    };

    //using Terminal = Array<Point, 2>; 
    using Terminal = std::vector<Point>;

    class Segment
    {
        public:

        Segment() = default;
        Segment(const Point&, const Point&);
        Segment(const Segment& other);

        vec3<double> normal() const;
        //vec3<double> centroid() const;
        const Terminal& vertex() const;
        const Point& vertex(int i) const;
        //bool do_intersect(const Point&) const;
        //bool do_intersect(const Segment& s, bool point_on_segment_means_inter) const;
        void rotate_points(double ang, double axis, const vec3<double>& rot_axis);
        const vec3<double>& min() const;
        const vec3<double>& max() const;
        double min(int i) const;
        double max(int i) const;
        void move_points(const vec3<double>& v);
        double len() const;
        //vec3<double> area() const;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & vertex_;
            ar & min_;
            ar & max_;
        }

        private:
        
        Terminal vertex_;
        vec3<double> min_;
        vec3<double> max_;

        void set_bbox();

        friend class boost::serialization::access;
    };

    class Polygon
    {
        using Vertices = Array<Point, 4>; 
        using Edges = Array<Segment, 4>; 

        public:

        Polygon() = default;
        Polygon(const std::vector<Point>& vertex);

        bool degenerate() const;
        void set_vertices(const std::vector<Point>& vertex);
        void set_edge();
        void rotate_points(double ang, double axis, const vec3<double>& rot_axis);
        void move_points(const vec3<double>& v);
        const Vertices& vertex() const;
        const Point& vertex(int i) const;
        const Edges& edge() const;
        const Segment& edge(int i) const;
        double signed_area() const;
        vec3<double> centroid() const;
        bool do_intersect(const Segment&) const;
        //bool do_intersect(const Point&) const;
        //bool do_intersect(const vec3<double>& r) const;
        //bool do_intersect(const Polygon&) const;
        double longest_edge_length() const;
        vec3<double> normal() const;
        double normal(int i) const;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & vertex_;
            ar & edge_;
            ar & min_;
            ar & max_;
            //ar & nvertex_;
        }

        private:

        vec3<double> min_;
        vec3<double> max_;
        Vertices vertex_;
        Edges edge_;
        mutable boost::optional<double> signed_area_;
        //int nvertex_;
        mutable boost::optional<vec3<double>> normal_; // don't use this directly even for internal usage. use normal().
        void set_bbox();
        std::vector<double> centroid_ab(int i0, int i1) const;

        friend class boost::serialization::access;
    };

    /*class Face
    {
        using Vertices = Array<Point, 8>; 

        public:

        void set_vertices(const std::vector<Point>& vertex);

        private:

        Vertices vertex_;
    };*/

    enum Shape
    {
        tet,
        hex,
        pri,
        quad,
        tri,
        undef,
    };

    class Polyhedron
    {
        using Vertices = Array<Point, 8>; 
        using Faces = Array<Polygon, 6>; 

        public:

        Polyhedron() = default;
        Polyhedron(const std::vector<Point>& vertex, Shape shape);

        bool degenerate() const;
        Shape shape() const;
        void rotate_points(double ang, double axis, const vec3<double>& rot_axis);
        void move_points(const vec3<double>& v);
        const Vertices& vertices() const;
        const Faces& faces() const;
        vec3<double> centroid() const;
        bool do_intersect(const vec3<double>&) const;
        bool do_intersect(const Segment&) const;
        double volume() const;
        const vec3<double>& min() const;
        const vec3<double>& max() const;
        double min(int i) const;
        double max(int i) const;
        const Polygon& face(int i) const;
        void set_min(vec3<double> min);
        void set_max(vec3<double> max);
        void set_bbox(vec3<double> min, vec3<double> max);
        void set_bbox(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax);
        void set_shape(Shape shape);
        void set_faces();
        void set_vertices_from_bbox();
        size_t mem() const;

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & vertex_;
            ar & face_;
            ar & min_;
            ar & max_;
            ar & shape_;
        }

        private:

        Vertices vertex_;
        Faces face_;
        mutable boost::optional<double> volume_;
        Shape shape_;
        //int nvertex_;
        vec3<double> min_; 
        vec3<double> max_; 

        //double volume_triface(std::array<Polygon>::const_iterator begin, std::array<Polygon>::const_iterator end) const;
        double volume_triface(const Polygon* begin, const Polygon* end) const;
        void set_bbox();
        void set_vertices(const std::vector<Point>& vertex);
        //vec3<double> centroid(std::array<Polygon>::const_iterator begin, std::array<Polygon>::const_iterator end) const;

        friend class boost::serialization::access;
    };

    class RegularHexahedron: public Polyhedron
    {
        // Opposing faces are of equal size.
        // Dimensions can be different.

        public:

        RegularHexahedron();
        RegularHexahedron(vec3<double> min, vec3<double> max);
        RegularHexahedron(double minx, double miny, double maxx, double maxy, double zmin, double zmax);
        RegularHexahedron(const std::vector<Point>& points);
        RegularHexahedron(const Polyhedron&);

        double volume_of_overlap(const RegularHexahedron& a);
        bool extend(const RegularHexahedron& other);
        using Polyhedron::do_intersect;
        bool do_intersect(const RegularHexahedron&) const;
        bool do_contain(const RegularHexahedron& other) const;

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Polyhedron>(*this);
        }

        private:

        void set_bbox_from_points(const std::vector<Point>& points);
        friend class boost::serialization::access;
    };

    using AABB = RegularHexahedron;

    class Outline
    {
        public:

            /*void build(const std::vector<std::vector<vec3<double>>>& point, int dummyrank);
            //bool contain(vec3<double> p, bool strict) const;
            const Tag& tag() const;
            void set_tag(const Tag& t);
            const std::vector<Segment>& segment() const;
            //bool do_intersect(const std::vector<vec3<double>>& vertices) const;
            bool do_contain(const vec3<double>& v, bool strict) const;
            //Polygon_2 polygon() const;
            void build_polygon(int dummyrank);*/

            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & tag_;
                ar & edge_;
                //ar & polygon_;
            }
            
        private:
            friend class boost::serialization::access;
            Tag tag_;
            std::vector<Segment> edge_;

    };

    int orientation (const vec3<double>& p, const vec3<double>& q, const vec3<double>& r);
    Polyhedron create_with_sweep(const Polygon& polygon, const vec3<double>& v);
    Polygon create_with_sweep(const Segment& segment, const vec3<double>& v);
}

#endif
