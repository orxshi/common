#include "boundary.h"

namespace Common
{
    std::map<int, std::string> BT::type_ = []
    {
        std::map<int, std::string> type_;
        BT::type_[-1] = "undefined";
        BT::type_[1] = "wall";
        BT::type_[2] = "dirichlet";
        BT::type_[3] = "empty";
        BT::type_[4] = "interior";
        BT::type_[9] = "farfield";
        BT::type_[10] = "partition";
        return type_;
    }();

    int BT::to_int(std::string s)
    {
        bool found = false;
        int i;
        for (auto it = type_.begin(); it != type_.end(); ++it)
        {
            if (it->second == s)
            {
                i = it->first;
                found = true;
                break;
            }
        }

        assert(found);
        return i;
    }

    std::string BT::to_string(int i)
    {
        auto it = type_.find(i);
        assert(it != type_.end());
        return it->second;
    }


    boundary_t::boundary_t()
    {
        type_ = std::make_pair(BT::to_int("undefined"), "undefined");
    }

    boundary_t::boundary_t(int i)
    {
        type_ = std::make_pair(i, BT::to_string(i));
    }

    boundary_t::boundary_t(std::string s)
    {
        type_ = std::make_pair(BT::to_int(s), s);
    }

    bool boundary_t::operator==(std::string s) const
    {
        return (type_.second == s);
    }

    bool boundary_t::operator!=(std::string s) const
    {
        return !(*this == s);
    }

    bool boundary_t::operator==(const boundary_t& other) const
    {
        return ((type_.first == other.type_.first) && (type_.second == other.type_.second));
    }

    bool boundary_t::operator!=(const boundary_t& other) const
    {
        return !(*this == other);
    }

    std::ostream& operator<<(std::ostream& os, const boundary_t& bt)
    {
        os << bt.type_.second;
        return os;
    }
}
