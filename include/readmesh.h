#ifndef COMMON_READMESH_H
#define COMMON_READMESH_H

#include "commonmesh.h"

namespace Common
{
    enum class gmsh
    {
        tri = 2,
        quad = 3,
        tet = 4,
        hex = 5,
        pri = 6,
    };

    void read_mesh(Mesh& farfield, Mesh& interior, Mesh& wall, Mesh& dirichlet, Mesh& empty, const boost::mpi::communicator& comm, std::string file_name);
}

#endif
