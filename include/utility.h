#ifndef COMMON_UTILITY_H
#define	COMMON_UTILITY_H

#include <vector>
#include <algorithm>

namespace Common
{
    template<class T> void remove_merge_dup(std::vector<T>& vec)
    {
        std::sort(vec.begin(), vec.end());

        auto beg = vec.begin();
        while(true)
        {
            beg = std::adjacent_find(beg, vec.end());
            auto temp_address = &(*beg);
            if (beg == vec.end()) return;

            beg->merge(*std::next(beg));
            vec.erase(std::next(beg));
            assert(&(*beg) == temp_address);
        }

        auto it = std::adjacent_find(vec.begin(), vec.end());
        assert(it == vec.end());
    }

    template<class T> void remove_dup(std::vector<T>& vec)
    {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }
}

#endif
