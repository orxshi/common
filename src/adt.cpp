#include "adt.h"

namespace Common
{
    /*bool ADT::extend_aabb(const AABB& other_aabb)
    {
        AABB old_aabb = aabb_;
        bool b = aabb_.extend(other_aabb);

        assert(aabb_.min(0) <= old_aabb.min(0));
        assert(aabb_.min(1) <= old_aabb.min(1));
        assert(aabb_.max(0) >= old_aabb.max(0));
        assert(aabb_.max(1) >= old_aabb.max(1));

        assert(old_aabb.min(0) != old_aabb.max(0));
        assert(old_aabb.min(1) != old_aabb.max(1));

        assert(aabb_.min(0) != aabb_.max(0));
        assert(aabb_.min(1) != aabb_.max(1));

        return b;
    }*/

    void ADT::rotate(double ang, int axis, const vec3<double>& rot_point)
    {
        assert(root_ != NULL);

        root_->rotate(ang, axis, rot_point);
    }

    void ADT::move(const vec3<double>& v)
    {
        assert(root_ != NULL);

        root_->move(v);
    }

    ADT::ADT(const ADT& other): root_(NULL)
    {
        //std::cout << "copy cons" << std::endl;
        if (other.root_ != NULL)
        {
            //delete root_;
            root_ = new Node(*other.root_);

            root_->construct_id_address_map(id_address_map);
        }
        aabb_ = other.aabb();
        assert(!aabb_.faces().empty());
        //std::cout << "copied cons" << std::endl;
    }

    ADT& ADT::operator=(const ADT& other)
    {
        destroy_tree();

        if (other.root_ != NULL)
        {
            root_ = new Node(*other.root_);
            root_->construct_id_address_map(id_address_map);
        }
        aabb_ = other.aabb();
        return *this;
    }

    ADT::ADT(ADT&& other)
    {
        //std::cout << "move cons" << std::endl;
        root_ = other.root_;
        other.root_ = NULL;

        id_address_map = other.id_address_map;
        other.id_address_map.clear();
        aabb_ = other.aabb();
        assert(!aabb_.faces().empty());
    }

    /*ADT::ADT(const std::vector<ADTPoint>& points): ADT()
      {
      init_root(points);
      insert(points);
      }*/

    ADT::ADT(const std::vector<ADTPoint>& points): root_(NULL)
    {
        init_root(points);
        assert(!aabb_.faces().empty());
        root_->add_point(points.back());
        id_address_map.insert(std::make_pair(points.back().idx(), root_));
        insert(points);
        assert(root_->p() != NULL);
    }

    ADT::ADT(): root_(NULL) {}

    ADT::ADT(const std::vector<double>& dim): root_(NULL)
    {
        init_root(dim);
        assert(!aabb_.faces().empty());
    }

    ADT::ADT(const std::vector<double>& dim, const std::vector<ADTPoint>& points): root_(NULL)
    {
        init_root(dim);
        assert(!aabb_.faces().empty());
        root_->add_point(points.back());
        id_address_map.insert(std::make_pair(points.back().idx(), root_));
        insert(points);
        assert(root_->p() != NULL);
    }

    ADT::~ADT()
    {    
        destroy_tree();
    }

    void ADT::destroy_tree()
    {   
        //std::cout << "destroying tree" << std::endl;
        //std::cout << "root_ = " << root_ << std::endl;
        delete root_;
        root_ = NULL;
        //std::cout << "destroyed tree" << std::endl;
    }

    void ADT::init_root(const std::vector<ADTPoint>& points)
    {
        assert(root_ == NULL);

        std::vector<double> c(ADT_VAR, BIG_POS_NUM);
        std::vector<double> d(ADT_VAR, BIG_NEG_NUM);

        for (const ADTPoint& p: points)
        {
            for (int i=0; i<ADT_VAR; ++i)
            {
                c[i] = std::min(p.dim(i), c[i]);
                d[i] = std::max(p.dim(i), d[i]);
            }
        }

        root_ = new Node(0, c, d);
        aabb_.set_bbox(c[0], c[2], c[4], d[1], d[3], d[5]);
        aabb_.set_vertices_from_bbox();
        aabb_.set_faces();
        assert(aabb_.min(0) == root_->c(0));
        assert(aabb_.min(1) == root_->c(2));
        assert(aabb_.max(0) == root_->d(1));
        assert(aabb_.max(1) == root_->d(3));

        assert(aabb_.min(0) != aabb_.max(0));
        assert(aabb_.min(1) != aabb_.max(1));
        //root_ = new Node(0, c, d, points.back());
        //id_address_map.insert(std::make_pair(points.back().idx(), root_));
    }

    const AABB& ADT::aabb() const
    {
        return aabb_;
    }

    void ADT::init_root(std::vector<double> dim)
    {
        assert(root_ == NULL);
        assert(dim.size() == ADT_VAR);

        //dim.push_back(0.); // zmin.
        //dim.push_back(0.); // zmax.

        root_ = new Node(0, dim, dim);
        aabb_.set_bbox(dim[0], dim[2], dim[4], dim[1], dim[3], dim[5]);
        aabb_.set_vertices_from_bbox();
        aabb_.set_faces();
        assert(aabb_.min(0) == root_->c(0));
        assert(aabb_.min(1) == root_->c(2));
        assert(aabb_.max(0) == root_->d(1));
        assert(aabb_.max(1) == root_->d(3));

        assert(aabb_.min(0) != aabb_.max(0));
        assert(aabb_.min(1) != aabb_.max(1));
    }

    bool ADT::insert(const ADTPoint& point)
    {
        assert(root_ != NULL);
        //assert(point.dim(0) < 601); // delete laster;
        //assert(point.dim(1) < 601); // delete laster;
        //assert(point.dim(2) < 601); // delete laster;
        //assert(point.dim(3) < 601); // delete laster;

        int count = id_address_map.count(point.idx());
        if (count != 0)
            return false;

        if (root_->p() == NULL)
        {
            root_->add_point(point);
            id_address_map.insert(std::make_pair(point.idx(), root_));
            assert(root_->p() != NULL);
            return true;
        }
        else
        {
            Node* node = root_->insert(point);
            if (node != NULL)
            {
                id_address_map.insert(std::make_pair(point.idx(), node));
                return true;
            }
        }

        return false;
    }

    const Node* const ADT::root() const
    {
        return root_;
    }

    void ADT::search(const Node* node, const ADTPoint& targetPoint, std::vector<const Node*>& searchStack, std::vector<int>& id, bool verbose) const
    {
        if (node != NULL)
        {
            //assert(node->p() != NULL); // not sure about this.
            if (verbose)
            {
                if (node->p() != NULL)
                    std::cout << "idx = " << node->p()->idx() << std::endl;
            }
            if (node->overlap(targetPoint, verbose))
            {
                if (verbose)
                    std::cout << "node overlaps target" << std::endl;
                id.push_back(node->p()->idx());

                node->searchChildren(targetPoint, searchStack);

                if (!searchStack.empty())
                {
                    const Node* last = searchStack.back();   
                    //searchStack.back() = NULL;                
                    searchStack.pop_back();
                    search(last, targetPoint, searchStack, id);
                }
            }
            else
            {
                if (verbose)
                    std::cout << "searching children" << std::endl;
                node->searchChildren(targetPoint, searchStack);

                if (!searchStack.empty())
                {
                    if (verbose)
                        std::cout << "search stack is NOT empty" << std::endl;
                    const Node* last = searchStack.back();        
                    searchStack.pop_back();
                    search(last, targetPoint, searchStack, id, verbose);
                }
                else
                {
                    if (verbose)
                        std::cout << "search stack is empty" << std::endl;
                }
            }
        }
        else
        {
            if (verbose)
                std::cout << "node is null" << std::endl;
        }
    }

    std::vector<int> ADT::search(const ADTPoint& targetPoint, bool verbose) const
    {
        //bool verbose = true;

        //if (targetPoint.idx() == 522)
            //verbose = true;
        //else
            //verbose = false;

        if (verbose)
            std::cout << "starting adt search" << std::endl;

        std::vector<const Node*> searchStack;    
        std::vector<int> id;    

        assert(root_ != NULL);
        //bool overlap = root_->doRegionOverlapCheckAll(targetPoint);
        Point p(targetPoint.dim(0), targetPoint.dim(2), targetPoint.dim(4));
        //bool overlap = aabb_.do_intersect(p);
        assert(!aabb_.faces().empty());
        bool overlap = aabb_.do_intersect(p.r());
        if (overlap)
        {
            if (verbose)
            {
                std::cout << "adt region overlaps" << std::endl;
                //std::cout << "has it = " << query(660) << std::endl;
            }
            //assert(root_->p() != NULL); // root.p might be NULL after disconneting root.p from sp.
            //std::cout << "root idx = " << root_->p()->idx() << std::endl;
            //search(root_, targetPoint, searchStack, id, verbose);
            search(root_, targetPoint, searchStack, id, verbose);
        }
        else
        {
            if (verbose)
            {
                std::cout << "adt region does not overlap" << std::endl;
                std::cout << "c0 = " << root_->c(0) << std::endl;
                std::cout << "c2 = " << root_->c(2) << std::endl;
                std::cout << "c4 = " << root_->c(4) << std::endl;
                std::cout << "d1 = " << root_->d(1) << std::endl;
                std::cout << "d3 = " << root_->d(3) << std::endl;
                std::cout << "d5 = " << root_->d(5) << std::endl;
                std::cout << "dim = " << targetPoint.dim(0) << std::endl;
                std::cout << "dim = " << targetPoint.dim(2) << std::endl;
                std::cout << "dim = " << targetPoint.dim(4) << std::endl;
                std::cout << "aabb = " << aabb_.min(0) << std::endl;
                std::cout << "aabb = " << aabb_.min(1) << std::endl;
                std::cout << "aabb = " << aabb_.min(2) << std::endl;
                std::cout << "aabb = " << aabb_.max(0) << std::endl;
                std::cout << "aabb = " << aabb_.max(1) << std::endl;
                std::cout << "aabb = " << aabb_.max(2) << std::endl;
                std::cout << "p = " << p.r(0) << std::endl;
                std::cout << "p = " << p.r(1) << std::endl;
                std::cout << "p = " << p.r(2) << std::endl;
            }
        }

        return id;
    }


    bool ADT::remove(int id)
    {
        auto it = id_address_map.find(id);
        if (it == id_address_map.end()) return false;

        Node* node = it->second;

        if (node == NULL) return false;

        if(node->p() == NULL)
        {
            std::cout << "id = " << id << std::endl;
        }
        assert(node->p() != NULL);

        node->remove_point();

        id_address_map.erase(it);
        auto itt = id_address_map.find(id);
        assert(itt == id_address_map.end());

        return true;
    }

    bool ADT::query(int id) const
    {
        auto it = id_address_map.find(id);
        if (it == id_address_map.end()) return false;

        assert(it->second->p() != NULL);
        assert(it->second->p()->idx() == id);

        return true;
    }

    void ADT::insert(const std::vector<ADTPoint>& points)
    {
        assert(!points.empty());


        for (auto p=points.begin(); p!=points.end(); ++p)
        {
            bool success = insert(*p);
            assert(success || p == points.end() - 1);
        }
    }
}
