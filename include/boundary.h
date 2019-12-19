#ifndef COMMON_BOUNDARY_H
#define COMMON_BOUNDARY_H

#include <map>
#include <cassert>

namespace Common
{
    struct BT
    {
        static std::map<int, std::string> type_;

        //BT();

        static int to_int(std::string s);
        static std::string to_string(int i);
    };

    struct boundary_t
    {
        std::pair<int, std::string> type_;

        boundary_t();
        boundary_t(int i);
        boundary_t(std::string s);
        bool operator==(std::string s) const;
        bool operator!=(std::string s) const;
        bool operator==(const boundary_t& other) const;
        bool operator!=(const boundary_t& other) const;

        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & type_;
        }
    };

    std::ostream& operator<<(std::ostream& os, const boundary_t& bt);
}

#endif
