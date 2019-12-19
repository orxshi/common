#include "sp.h"

namespace Common
{
    void SpatialPartition::connect_cells_uniproc()
    {
        for (Mesh& m: mesh_)
        {
            m.connect_cells();
        }
    }

    void SpatialPartition::connect_cells()
    {
        if (!master())
        {
            for (Mesh& m: mesh_)
            {
                m.connect_cells();
            }
        }
    }
}

