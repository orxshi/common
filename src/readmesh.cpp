#include "readmesh.h"

namespace Common
{
    std::ifstream& go_to_beg_of_line(std::ifstream& file, int num)
    {   
        file.seekg(std::ios::beg);
        for(int i=0; i<num-1; ++i)
        {
            file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        }

        return file;
    }

    void read_mesh_GMSH(Mesh& mesh, std::string file_name)
    {
        std::string temps;
        int tempi;
        int line_number = 1;
        int tag, igs, n_tags, geo, n_part, n_total, n_point, part_tag;
        gmsh gs;
        boundary_t phys("undefined");

        std::ifstream in;
        in.open (file_name);
        if (!in.is_open())
        {
            std::cout << "file " << file_name << " could not be opened." << std::endl;
            return;
        }
        assert(in.is_open());

        in >> temps; // mesh format.
        in >> temps; // version.
        assert(temps == "2.2"); // msh2 format.
        in >> temps; // version.
        in >> temps; // version.
        in >> temps; // end mesh format.
        in >> temps; // nodes.
        in >> n_point; // number under "$Nodes" is the number of points.
        line_number += 5;

        //point_.reserve(n_point);
        for (int i=0; i<n_point; ++i)
        {
            int ptag;
            double x, y, z;
            in >> ptag;
            in >> x;
            in >> y;
            in >> z;
            ++line_number;
            //assert(ptag != 1449);

            //MeshPoint mp;
            //mp.set_tag(ptag);
            //mp.set_parent_mesh(tag_);
            //point_.push_back(std::move(mp));
            mesh.add_point(MeshPoint(x, y, z), ptag, n_point);
        }

        //mesh.set_point_tag_index_map();
        //mesh.shrink_points(); // uncomment
        mesh.sort_points(); // uncomment

        in >> temps; // end of nodes.
        in >> temps; // elements.
        in >> n_total; // the number under "$Elements" is total number of elements which includes boundary faces and cells.
        line_number += 3;

        // read elements.
        int n_bface = 0;
        for (int e=0; e<n_total; ++e)
        {
            in >> tag; // tag of boundary face or cell.
            in >> igs; // geometric shape of boundary face or cell.
            gs = static_cast<gmsh>(igs);
            in >> n_tags; // number of GMSH tags.

            int tag_count = 0;

            // read GMSH tags.
            while (n_tags > 0)
            {
                // read physical number.
                in >> tempi;
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> geo; // geometrical number.
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> n_part; // number of partitions to which element belongs.
                for (int i=0; i<n_part; ++i)
                {
                    in >> temps;
                }

                break;
            }

            if (gs == gmsh::tri)
            {
                in >> tempi;
                in >> tempi;
                in >> tempi;
            }
            else if (gs == gmsh::quad)
            {
                in >> tempi;
                in >> tempi;
                in >> tempi;
                in >> tempi;
            }
            else
            {
                break;
            }

            ++n_bface;
        }

        go_to_beg_of_line(in, line_number);

        for (int e=0; e<n_total; ++e)
        {
            in >> tag; // read tag of boundary face or cell.
            in >> igs; // read geometric shape of boundary face or cell.
            gs = static_cast<gmsh>(igs);
            in >> n_tags; // read number of GMSH tags.

            int tag_count = 0;

            while (n_tags > 0)
            {
                // read physical number.
                in >> tempi;
                phys = boundary_t(tempi);
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> geo; // read geometrical number.
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> n_part; // read number of partitions to which element belongs.
                //assert(n_part == 1);
                for (int i=0; i<n_part; ++i)
                {
                    in >> part_tag;
                }

                break;
            }

            if (e < n_bface)
            {
                Shape shape = Shape::undef;
                int n_vertex;
                switch (gs)
                {
                    case gmsh::tri:
                        n_vertex = 3;
                        shape = Shape::tri;
                        break;
                    case gmsh::quad:
                        n_vertex = 4;
                        shape = Shape::quad;
                        break;
                    default:
                        std::cout << "invalid n_vertexx" << std::endl;
                        std::cout << "gs = " << igs << std::endl;
                        std::cout << "tag = " << tag << std::endl;
                        std::cout << "n_tags = " << n_tags << std::endl;
                        exit(0);
                }
                std::vector<int> vtx;
                vtx.reserve(n_vertex);
                for (int i=0; i<n_vertex; ++i)
                {
                    in >> tempi;
                    vtx.push_back(tempi);
                    if (i != 0)
                    {
                        assert(vtx[i] != vtx[i-1]);
                    }
                }
                std::vector<MeshPoint> mpts;
                for (int z=0; z<vtx.size(); ++z)
                {
                    mpts.push_back(mesh.point(Tag(vtx[z])));
                }

                mesh.add_cell_only(MeshCell(Tag(tag), mesh.tag(), mpts, phys, shape));

                for (int i=0; i<mesh.cell().back().point().size(); ++i)
                {
                    if (i != 0)
                    {
                        assert(mesh.cell().back().point(i).tag() != mesh.cell().back().point(i-1).tag());
                    }
                }
            }
            else
            {
                Shape shape;
                int n_vertex;
                switch (gs)
                {
                    case gmsh::tet:
                        n_vertex = 4;
                        shape = Shape::tet;
                        break;
                    case gmsh::hex:
                        n_vertex = 8;
                        shape = Shape::hex;
                        break;
                    case gmsh::pri:
                        n_vertex = 6;
                        shape = Shape::pri;
                        break;
                    default:
                        std::cout << "igs: " << igs << std::endl;
                        std::cout << "invalid n_vertex" << std::endl;
                        std::cout << "igs = " << igs << std::endl;
                        std::cout << "tag = " << tag << std::endl;
                        std::cout << "n_tags = " << n_tags << std::endl;
                        exit(0);
                }

                // vertices.
                std::vector<int> vtx;
                vtx.reserve(n_vertex);

                for (int i=0; i<n_vertex; ++i)
                {
                    in >> tempi;
                    vtx.push_back(tempi);
                    if (i != 0)
                    {
                        assert(vtx[i] != vtx[i-1]);
                    }
                }

                std::vector<MeshPoint> mpts;

                for (int z=0; z<vtx.size(); ++z)
                {
                    mpts.push_back(mesh.point(Tag(vtx[z])));
                }

                mesh.add_cell_only(MeshCell(Tag(tag), mesh.tag(), mpts, phys, shape), n_total);
                for (int i=0; i<mesh.cell().back().point().size(); ++i)
                {
                    if (i != 0)
                    {
                        assert(mesh.cell().back().point(i).tag() != mesh.cell().back().point(i-1).tag());
                    }
                }

            }
        }

        in.close();

        mesh.sort_cells();
        mesh.update_points_from_cell_vertices();
    }

    void read_wall(Mesh& mesh, std::string file_name)
    {
        file_name.append("_wall.msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_dirichlet(Mesh& mesh, std::string file_name)
    {
        file_name.append("_dirichlet.msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_empty(Mesh& mesh, std::string file_name)
    {
        file_name.append("_empty.msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_farfield(Mesh& mesh, std::string file_name)
    {
        file_name.append("_farfield.msh");
        read_mesh_GMSH(mesh, file_name);
    }

    void read_interior_cells(Mesh& mesh, std::string file_name, int rank, bool uniproc)
    {
        if (!uniproc)
        {
            file_name.append("_");
            file_name.append(std::to_string(rank));
        }
        file_name.append(".msh");
        std::cout << "reading " << file_name << std::endl;
        read_mesh_GMSH(mesh, file_name);
    }

    void read_mesh(Mesh& farfield, Mesh& interior, Mesh& wall, Mesh& dirichlet, Mesh& empty, const boost::mpi::communicator& comm, std::string file_name)
    {
        bool uniproc = comm.size() == 1 ? true : false;
        int rank = comm.rank();
        assert(file_name != "");

        if (uniproc || rank != 0) {
            read_wall(wall, file_name);
        }

        std::cout << "read wall" << std::endl;

        if (uniproc || rank != 0) {
            read_dirichlet(dirichlet, file_name);
        }
        std::cout << "read diri" << std::endl;

        if (uniproc || rank != 0) {
            read_empty(empty, file_name);
        }
        std::cout << "read empty" << std::endl;

        if (uniproc || rank != 0) {
            read_farfield(farfield, file_name);
        }
        std::cout << "read farfield" << std::endl;

        if (uniproc || rank != 0) {
            read_interior_cells(interior, file_name, rank, uniproc);
            interior.set_all_cells_as_interior();
        }
        std::cout << "read interior" << std::endl;

        if (uniproc || rank != 0)
        {
            interior.connect_add_bou_to_interior(wall, boundary_t("wall"));
            std::cout << "connect wall" << std::endl;
            interior.connect_add_bou_to_interior(farfield, boundary_t("farfield"));
            std::cout << "connect farfield" << std::endl;
            interior.connect_add_bou_to_interior(dirichlet, boundary_t("dirichlet"));
            std::cout << "connect diri" << std::endl;
            interior.connect_add_bou_to_interior(empty, boundary_t("empty"));
            std::cout << "connect empty" << std::endl;
        }

        if (uniproc || rank != 0)
        {
            assert(!interior.cell().empty());
        }
    }
}
