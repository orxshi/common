#include "adt.h"

namespace Common
{
    void Node::rotate(double ang, int axis, const vec3<double>& rot_point)
    {
        if (p_ != NULL)
            p_->rotate(ang, axis, rot_point);

        if (left_ != NULL)
            left_->rotate(ang, axis, rot_point);

        if (right_ != NULL)
            right_->rotate(ang, axis, rot_point);
    }

    void Node::move(const vec3<double>& v)
    {
        if (p_ != NULL)
            p_->move(v);

        if (left_ != NULL)
            left_->move(v);

        if (right_ != NULL)
            right_->move(v);
    }

    double Node::c(int i) const
    {
        assert(i >= 0);
        assert(i < ADT_VAR);
        return c_[i];
    }

    double Node::d(int i) const
    {
        assert(i >= 0);
        assert(i < ADT_VAR);
        return d_[i];
    }

    void Node::construct_id_address_map(std::map<int, Node*>& map)
    {
        if (p_ != NULL)
        {
            map.insert(std::make_pair(p_->idx(), this));
        }

        if (left_ != NULL)
            left_->construct_id_address_map(map);
        if (right_ != NULL)
            right_->construct_id_address_map(map);
    }

    Node::Node(const Node& other): c_(other.c_), level_(other.level_), d_(other.d_), key_(other.key_), left_(NULL), dim_(other.dim_), p_(NULL), right_(NULL)
    {
        //std::cout << "copy cons of node" << std::endl;
        if (other.p_ != NULL)
        {
            //delete p_;
            p_ = new ADTPoint(*other.p_);
        }

        if (other.left_ != NULL)
        {
            //delete left_;
            left_ = new Node(*other.left_);
        }

        if (other.right_ != NULL)
        {
            //delete right_;
            right_ = new Node(*other.right_);
        }
        //std::cout << "exit copy cons of node" << std::endl;
    }

    Node::Node(unsigned int level, const std::vector<double>& c, const std::vector<double>& d): left_(NULL), p_(NULL), right_(NULL), c_(c), level_(level), d_(d)
    {
        dim_ = level_ % ADT_VAR;
        key_ = 0.5 * (c[dim_] + d[dim_]);
    }

    Node::Node(unsigned int level, const std::vector<double>& c, const std::vector<double>& d, const ADTPoint& p): Node(level, c, d)
    {
        p_ = new ADTPoint(p);
    }

    Node::~Node()
    {    
        destroy_children();
    }

    void Node::destroy_children()
    {
        if (p_ != NULL)
            delete p_;

        if (left_ != NULL)
            left_->destroy_children();

        if (right_ != NULL)
            right_->destroy_children();
    }

    unsigned int Node::level() const
    {
        return level_;
    }

    double Node::key() const
    {
        return key_;
    }

    Node* Node::insert(const ADTPoint& point)
    {
        Node* node = NULL;

        if (point.dim(dim_) < key_)
        {
            bool added = insert_left(point);
            if (added)
                return left_;
            node = left_->insert(point);
        }
        else
        {
            bool added = insert_right(point);
            if (added)
                return right_;
            node = right_->insert(point);
        }

        return node;
    }

    void Node::add_point(const ADTPoint& point)
    {
        p_ = new ADTPoint(point);
        //isEmpty_ = false;
    }

    bool Node::insert_left(const ADTPoint& point)
    {
        if (left_ == NULL)
        {
            std::vector<double> newd = d_;

            for (int d=0; d<ADT_VAR; ++d)
            {
                if (d == dim_)
                {
                    newd[d] = key_;
                }
            }

            left_ = new Node(level_ + 1, c_, newd);
        }

        if (left_->p_ == NULL)
        {
            left_->add_point(point);
            return true;
        }

        return false;
    }

    bool Node::insert_right(const ADTPoint& point)
    {
        if (right_ == NULL)
        {
            std::vector<double> newc = c_;

            for (int d=0; d<ADT_VAR; ++d)
            {
                if (d == dim_)
                {
                    newc[d] = key_;
                }
            }

            right_ = new Node(level_ + 1, newc, d_);
        }

        if (right_->p_ == NULL)
        {
            right_->add_point(point);
            return true;
        }

        return false;
    }

    bool Node::doRegionOverlapCheckAll(const ADTPoint& targetPoint) const
    {
        assert(false);
        /*//std::cout << "checck all" << std::endl;
        assert(c_.size() == ADT_VAR);
        assert(d_.size() == ADT_VAR);
        assert(targetPoint.dim().size() == ADT_VAR);

        for (int d=0; d<TAILOR_ADT_DIM; ++d)
        {
            double dis1 = c_[2*d] - targetPoint.dim(2*d+1);
            double dis2 = d_[2*d+1] - targetPoint.dim(2*d);
            //if (c_[2*d] > targetPoint.dim(2*d+1))
            if (dis1 > ZERO)
                return false;
            //if (d_[2*d+1] < targetPoint.dim(2*d))
            if (dis2 < ZERO)
                return false;
        }
        //std::cout << "end checck all" << std::endl;

        return true;*/

        //Point p(targetPoint.dim(0), targetPoint.dim(2));
        //return aabb_.do_intersect(p);
    }

    bool Node::doRegionOverlap(const ADTPoint& targetPoint) const
    {
        if ((dim_ % 2) == 0)
        {
            double dis1 = c_[dim_] - targetPoint.dim(dim_+1);
            //if (c_[dim_] > targetPoint.dim(dim_+1))
            if (dis1 > ZERO)
                return false;
        }
        else
        {
            double dis1 = d_[dim_] - targetPoint.dim(dim_-1);
            //if (d_[dim_] < targetPoint.dim(dim_-1))
            if (dis1 < ZERO)
                return false;
        }

        return true;
    }

    bool Node::overlap(const ADTPoint& targetPoint, bool verbose) const
    {
        //bool verbose = false;

        //if (targetPoint.idx() == 537)
        //verbose = true;
        //else
        //verbose = false;

        //if (p_ == NULL) return false;
        //if (!doCubesOverlap(targetPoint)) return false;
        //if (!compareFunction(targetPoint)) return false;

        if (p_ != NULL)
        {
            if (verbose)
                std::cout << "p is not null" << std::endl;
            if (p_->overlap(targetPoint, verbose)) return true;
        }
        else
        {
            if (verbose)
                std::cout << "p is null" << std::endl;
        }

        return false;
    }

    void Node::searchChildren(const ADTPoint& targetPoint, std::vector<const Node*>& searchStack) const
    {
        if (left_ != NULL)
        {
            if (left_->doRegionOverlap(targetPoint))
            {
                searchStack.push_back(left_);
            }
        }

        if (right_ != NULL)
        {
            if (right_->doRegionOverlap(targetPoint))
            {
                searchStack.push_back(right_);
            }
        }
    }

    const ADTPoint* Node::p() const
    {
        return p_;
    }

    void Node::remove_point()
    {
        assert(p_ != NULL);
        delete p_;
        p_ = NULL;
    }
}
