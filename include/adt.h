#ifndef COMMON_ADT_H
#define	COMMON_ADT_H

#include <map>
#include "vec3.h"
#include "commoncore.h"
#include "geom.h"
#include "rotationmatrix.h"
//#include "aabb.h"

//#define COMMON_ADT_DIM NDIM
//#define COMMON_ADT_VAR (2*COMMON_ADT_DIM)

namespace Common
{
    const int ADT_DIM = NDIM;
    const int ADT_VAR = 2 * ADT_DIM;

    class ADTPoint
    {
        int idx_;
        std::vector<double> dim_;
        /*
         * 0: xmin
         * 1: xmax
         * 2: ymin
         * 3: ymax
         * 4: zmin
         * 5: zmax
         */

        public:

        ADTPoint();
        //ADTPoint(const std::vector<Point>& vertex, int idx);
        //ADTPoint(const std::vector<Point>& vertex): ADTPoint(vertex, -1) {}
        template<typename Rai> ADTPoint(Rai begin, Rai end, int idx);
        template<typename Rai> ADTPoint(Rai begin, Rai end): ADTPoint(begin, end, -1) {}
        ADTPoint(const vec3<double>& p, int idx);
        ADTPoint(const vec3<double>& p): ADTPoint(p, -1) {}

        void rotate(double ang, int axis, const vec3<double>& rot_point);
        void move(const vec3<double>& v);
        void set_idx(int i);
        double dim(int i) const;
        const std::vector<double>& dim() const;
        int idx() const;
        bool overlap(const ADTPoint& other, bool verbose=false) const;
    };

    class Node
    {
        std::vector<double> c_, d_;
        ADTPoint* p_;
        double key_;
        Node* left_;
        Node* right_;
        unsigned int level_;
        unsigned int dim_;

        void destroy_children();
        bool insert_left(const ADTPoint& point);
        bool insert_right(const ADTPoint& point);

        public:

        Node(unsigned int level, const std::vector<double>& c, const std::vector<double>& d);
        Node(unsigned int level, const std::vector<double>& c, const std::vector<double>& d, const ADTPoint& p);
        Node(const Node& other);
        ~Node();

        void rotate(double ang, int axis, const vec3<double>& rot_point);
        void move(const vec3<double>& v);
        double c(int i) const;
        double d(int i) const;
        void add_point(const ADTPoint& point);
        const ADTPoint* p() const;
        const Node* left() const;
        const Node* right() const;
        unsigned int level() const;
        unsigned int dim() const;
        double key() const; 
        Node* insert(const ADTPoint& point);
        bool doRegionOverlap(const ADTPoint& targetPoint) const;
        bool doRegionOverlapCheckAll(const ADTPoint& targetPoint) const;
        void searchChildren(const ADTPoint& targetPoint, std::vector<const Node*>& searchStack) const;
        bool overlap(const ADTPoint& targetPoint, bool verbose=false) const;
        void remove_point();
        void construct_id_address_map(std::map<int, Node*>& map);
    };

    class ADT
    {
        Node *root_;
        std::map<int, Node*> id_address_map;
        AABB aabb_;

        void destroy_tree();
        void search(const Node* node, const ADTPoint& targetPoint, std::vector<const Node*>& searchStack, std::vector<int>& ids, bool verbose=false) const;
        void init_root(std::vector<double> dim);
        void init_root(const std::vector<ADTPoint>& points);

        public:

        ADT(); 
        ADT(const std::vector<double>& dim);
        ADT(const std::vector<double>& dim, const std::vector<ADTPoint>& points);
        ADT(const std::vector<ADTPoint>& points);
        ADT(ADT&& other);
        ADT(const ADT& other);
        ADT& operator=(const ADT& other);
        ~ADT();

        const AABB& aabb() const;
        //bool extend_aabb(const AABB& other_aabb);
        void rotate(double ang, int axis, const vec3<double>& rot_point);
        void move(const vec3<double>& v);
        const Node* const root() const;
        bool remove(int id);
        void insert(const std::vector<ADTPoint>& points);
        bool insert(const ADTPoint& point);
        std::vector<int> search(const ADTPoint& targetPoint, bool verbose=false) const;
        bool query(int id) const;
    };

}

#include "adt.hpp"

#endif
