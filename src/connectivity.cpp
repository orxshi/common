#include "commonmesh.h"
#include "commonmeshpoint.h"
#include "commonmeshcell.h"

namespace Common
{
    void Mesh::connect_cells()
    {
        for (MeshFace& mf: face_)
        {
            assert(mf.btype() != "undefined");
            if (mf.btype() == "partition")
            {
                assert(mf.parent_cell().size() == 1);
            }
        }

        for (MeshCell& mc: cell_)
        {
            for (MeshFace& mf: mc.face_p())
            {

                if (mf.btype() == "undefined")
                {
                    //std::cout << mf.parent_cell().size() << std::endl;
                    assert(mf.parent_cell().size() == 0);
                }
            }
        }






        auto addface = [&](MeshFace& mf, MeshFace& omf, int& facecounter)
        {
            MeshFace newface = mf;
            Tag newtag = Tag(facecounter);
            assert(newtag.isvalid());
            mf.set_tag(newtag);
            omf.set_tag(newtag);
            newface.set_tag(newtag);
            face_.push_back(newface);
            ++facecounter;
        };

        // this is needed after connection/disconnection of cells. but bad after connect_X_to_interior functions because there we set parent cells of faces. so call this during connection/dicconnection before connect cells. but anyway this is not needed for solver since no dis(connection) happens here.
        /*for (MeshCell& mc: cell_)
          {
          mc.remove_all_neighbors();
          for (MeshFace& mf: mc.face_p())
          {
          mf.remove_parent_cells();
          }
          }*/

        for (MeshCell& mc: cell_)
        {
            mc.set_boutype(boundary_t("undefined"));

            for (MeshFace& mf: mc.face_p())
            {
                if (mf.is_boundary())
                {
                    if (mf.parent_cell().size() != 2)
                    {
                        std::cout << "is_boundary: " << mf.is_boundary() << std::endl;
                        std::cout << "mf.parent_cell().size(): " << mf.parent_cell().size() << std::endl;
                        std::cout << "btype: " << mf.btype() << std::endl;
                    }
                    assert(mf.parent_cell().size() == 2);
                }
                if (mf.parent_cell().size() == 2) continue;
                assert(mf.parent_cell().size() != 1);

                if (std::count(mf.parent_cell().begin(), mf.parent_cell().end(), mc.tag()) == 0)
                {
                    mf.add_parent_cell(mc.tag());
                }
                assert(!mf.mesh_point().empty());
                const MeshPoint& p0 = point(mf.mesh_point(0));
                assert(!p0.parent_cell().empty());

                for (const Tag& nei: p0.parent_cell())
                {
                    assert(query(nei) != nullptr);
                    MeshCell& nc = cell_p(nei);
                    assert(nei.isvalid());
                    if (nei == mc.tag())
                    {
                        continue;
                    }

                    int counter = 0;
                    for (int i=0; i<mf.mesh_point().size(); ++i)
                    {
                        const MeshPoint& p = point(mf.mesh_point(i));
                        for (const Tag& _t: p.parent_cell())
                        {
                            if (_t == nei)
                            {
                                ++counter;
                                break;
                            }
                        }
                    }

                    if (counter == mf.mesh_point().size())
                    {
                        if (std::count(mc.pnei().begin(), mc.pnei().end(), nei) == 0)
                        {
                            mc.add_pnei(nei);
                        }
                        if (std::count(nc.pnei().begin(), nc.pnei().end(), mc.tag()) == 0)
                        {
                            nc.add_pnei(mc.tag());
                        }
                        mf.add_parent_cell(nei);
                        assert(mf.parent_cell().size() == 2);
                    }
                }

                //if (mf.parent_cell().size() != 2)
                //{
                    //std::cout << "is_boundary: " << mf.is_boundary() << std::endl;
                    //std::cout << "mf.parent_cell().size(): " << mf.parent_cell().size() << std::endl;
                //}
                //assert(mf.parent_cell().size() == 2);
            }
        }

        for (MeshCell& mc: cell_)
        {
            if (mc.near_boundary()) continue;
            if (mc.pnei().size() != mc.poly().faces().size())
            {
                mc.set_boutype(boundary_t("partition"));
            }
        }

        for (MeshCell& mc: cell_)
        {
            for (MeshFace& mf: mc.face_)
            {
                if (mf.parent_cell().size() == 1)
                {
                    assert(mf.btype() == "undefined");
                    mf.set_btype(boundary_t("partition"));
                }
            }
        }

        for (MeshCell& mc: cell_)
        {
            for (MeshFace& mf: mc.face_)
            {
                if (mf.parent_cell().size() == 2)
                {
                    if (mf.btype() == "undefined")
                    {
                        mf.set_btype(boundary_t("interior"));
                    }
                }
            }
        }

        assert(face_.empty());
        int facecounter = 0;

        for (MeshCell& mc: cell_)
        {
            for (MeshFace& mf: mc.face_)
            {
                if (mf.btype() != "interior")
                {
                    if (mf.parent_cell().size() == 2)
                    {
                        assert(!mf.tag().isvalid());
                        assert(mf.parent_cell().size() == 2);
                        const Tag& pc = mf.parent_cell(0);
                        MeshCell* nc = nullptr;
                        if (mf.btype() == "wall")
                        {
                            nc = &wall_boundary_p(pc);
                        }
                        else if (mf.btype() == "dirichlet")
                        {
                            nc = &dirichlet_boundary_p(pc);
                        }
                        else if (mf.btype() == "farfield")
                        {
                            nc = &farfield_boundary_p(pc);
                        }
                        else if (mf.btype() == "empty")
                        {
                            nc = &empty_boundary_p(pc);
                        }
                        else
                        {
                            std::cout << mf.btype() << std::endl;
                            assert(false);
                        }
                        MeshFace& omf = nc->face_[0];
                        addface(mf, omf, facecounter);
                    }
                    else if (mf.parent_cell().size() == 1)
                    {
                        const Tag& pc = mf.parent_cell(0);
                        MeshCell* nc = nullptr;
                        nc = &cell_p(pc);
                        MeshFace& omf = nc->face_[0];
                        addface(mf, omf, facecounter);
                    }
                    else
                    {
                        assert(false);
                    }
                }
                else
                {
                    if (mf.tag().isvalid())
                    {
                        continue;
                    }

                    assert(mf.parent_cell().size() == 2);
                    mf.set_btype(boundary_t("interior"));

                    bool neifound = false;
                    for (const Tag& pc: mf.parent_cell())
                    {
                        if (pc == mc.tag())
                        {
                            continue;
                        }

                        neifound = true;

                        MeshCell* nc = nullptr;
                        nc = &cell_p(pc);
                        assert(nc != nullptr);

                        bool facefound = false;
                        for (MeshFace& omf: nc->face_)
                        {
                            if (omf.parent_cell().size() != 2)
                            {
                                assert(omf.btype() == "partition");
                            }
                            for (const Tag& opc: omf.parent_cell())
                            {
                                if (opc == mc.tag())
                                {
                                    addface(mf, omf, facecounter);
                                    facefound = true;
                                    break;
                                }
                            }
                            if (facefound)
                            {
                                break;
                            }
                        }

                        assert(facefound);
                        if (!mf.tag().isvalid())
                        {
                            std::cout << "is boundary: " << mf.is_boundary() << std::endl;
                            std::cout << "btype: " << mf.btype() << std::endl;
                            std::cout << "nc face size: " << nc->face_.size() << std::endl;
                            std::cout << "nc face parent size: " << nc->face_[0].parent_cell().size() << std::endl;
                        }
                    }

                    if (!neifound)
                    {
                        std::cout << mf.parent_cell(0)() << std::endl;
                        std::cout << mf.parent_cell(1)() << std::endl;
                    }

                    assert(neifound);
                }

                assert(mf.tag().isvalid());
            }
        }

        for (const MeshCell& mc: cell_)
        {
            for (const MeshFace& mf: mc.face())
            {
                assert(mf.tag().isvalid());
                assert(!mf.parent_cell().empty());
            }
        }

        for (const MeshCell& mc: cell_)
        {
            assert(!mc.pnei().empty());
        }
    }
}
