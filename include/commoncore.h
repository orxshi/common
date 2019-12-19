#ifndef COMMON_CORE_H
#define	COMMON_CORE_H

namespace Common
{
    const double BIG_POS_NUM = 1e9;
    const double BIG_NEG_NUM = -1e9;
    const int MASTER = 0;
    const double ZERO = 1e-7;
    const int NDIM = 3;

    /*enum class boundary_t
    {
        undefined = -1,
        wall = 1,
        dirichlet = 2,
        empty = 3,
        interior = 4,
        periodicxl = 5,
        periodicxr = 6,
        periodicyl = 7,
        periodicyr = 8,
        farfield = 9,
        partition = 10,
    };*/

    enum OGA_cell_type_t
    {
        undefined = -1,
        non_resident = 0,
        receptor = 1,
        field = 2,
        hole = 3,
        mandat_receptor = 4,
        orphan = 5,
        hole_candidate = 6,
    };

    enum class LoadEstimType
    {
        area = 0,
        hybrid = 1,
        solver = 2,
        minmesh = 3,
    };
}

#endif
