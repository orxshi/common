#include "rm_index.h"

namespace Common
{
    RegularMeshIndex RegularMeshIndex::operator+(const RegularMeshIndex& other) const
    {
        return RegularMeshIndex(i_ + other.i(), j_ + other.j(), k_ + other.k());
    }

    void RegularMeshIndex::inci()
    {
        ++i_;
    }

    void RegularMeshIndex::incj()
    {
        ++j_;
    }

    void RegularMeshIndex::inck()
    {
        ++k_;
    }

    RegularMeshIndex::RegularMeshIndex(int i, int j, int k)
    {
        set(i, j, k);
    }

    void RegularMeshIndex::set(int i, int j, int k)
    {
        i_ = i;
        j_ = j;
        k_ = k;
    }

    const int RegularMeshIndex::i() const
    {
        return i_;
    }

    const int RegularMeshIndex::j() const
    {
        return j_;
    }

    const int RegularMeshIndex::k() const
    {
        return k_;
    }
}
