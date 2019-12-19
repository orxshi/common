#include "regular_mesh.h"
#include <iomanip>

namespace Common
{
    void RegularMesh::register_bincells(const std::vector<BinCell>& bincell, const std::deque<Mesh>& meshes)
    {
        size_t quarter = bincell.size() / 4;
        for (Bin& b: bin_)
        {
            b.reserve(quarter);
        }

        for (const BinCell& bc: bincell)
        {
            Tag mt = bc.mesh();
            Tag ct = bc.cell();

            AABB _aabb;
            const MeshCell* _cell = NULL;
            for (auto m=meshes.begin(); m!=meshes.end(); ++m)
            {
                if (mt == m->tag())
                {
                    assert(m->query(ct));
                    _cell = &m->cell(ct);
                    _aabb = AABB(_cell->poly());
                    break;
                }
            }

            assert(_cell != NULL);

            for (Bin& b: bin_)
            {
                if (b.aabb().do_intersect(_aabb))
                {
                    b.copy_bincell(bc);

                    //if (!_cell->is_ghost())
                    {
                        bool point_in_bin = b.aabb().do_intersect(_cell->poly().centroid());
                        if (point_in_bin)
                        {
                            b.increment_mesh_load(mt);
                        }
                    }
                }
            }
        }

        for (Bin& b: bin_)
        {
            b.shrink();
        }
    }

    void RegularMesh::register_resident_mesh(const Mesh& mesh, bool calcload, int rank)
    {
        //assert(rank != 0);
        // ghosts must be connected to interiors at this stage.
        assert(!mesh.cell().empty());
        for (const MeshCell& c: mesh.cell())
        {
            //if (mesh.tag()() == 0 && c.tag()() == 18003)
            //{
                //assert(false);
            //}
            assert(c.tag().isvalid());
            assert(mesh.tag().isvalid());

            assert(aabb_.do_intersect(AABB(c.poly())));

            {
                assert(aabb_.do_intersect(AABB(c.poly())));
                register_resident_cell(c, c.tag(), mesh.tag(), calcload, rank);
            }
        }
        {
            bool all_empty = true;
            for (const Bin& b0: bin())
            {
                if (!b0.cell().empty())
                {
                    all_empty = false;
                    break;
                }
            }
            assert(!all_empty);
        }
    }

    void RegularMesh::register_overlapping_mesh(const Mesh& mesh, bool calcload, int rank, bool adaptive)
    {
        assert(!mesh.cell().empty());
        for (const MeshCell& c: mesh.cell())
        {
            if (aabb_.do_intersect(AABB(c.poly())) == false) {
                continue;
            }
            assert(c.tag().isvalid());
            assert(mesh.tag().isvalid());

            assert(aabb_.do_intersect(AABB(c.poly())));

            for (const Bin& b: bin_)
            {
                assert(b.mesh_tag_index_map().size() == b.mesh_load().size());
            }

            {
                double dummy;
                register_cell(c, c.tag(), mesh.tag(), calcload, rank, adaptive, dummy, dummy, dummy, dummy, dummy);
            }

            for (const Bin& b: bin_)
            {
                assert(b.mesh_tag_index_map().size() == b.mesh_load().size());
            }
        }
        {
            bool all_empty = true;
            for (const Bin& b0: bin())
            {
                if (!b0.cell().empty())
                {
                    all_empty = false;
                    break;
                }
            }
            assert(!all_empty);
        }
    }

    void RegularMesh::register_mesh(const Mesh& mesh, bool calcload, int rank, bool adaptive, const std::shared_ptr<Profiler>& profiler)
    {
        assert(!mesh.cell().empty());

        double dur1 = 0.;
        double dur2 = 0.;
        double dur3 = 0.;
        double dur4 = 0.;
        double dur5 = 0.;
        for (const MeshCell& c: mesh.cell())
        {
            assert(c.tag().isvalid());
            assert(mesh.tag().isvalid());

            assert(aabb_.do_intersect(AABB(c.poly())));

            for (const Bin& b: bin_)
            {
                assert(b.mesh_tag_index_map().size() == b.mesh_load().size());
            }

            assert(aabb_.do_intersect(AABB(c.poly())));
            
            register_cell(c, c.tag(), mesh.tag(), calcload, rank, adaptive, dur1, dur2, dur3, dur4, dur5);

            for (const Bin& b: bin_)
            {
                assert(b.mesh_tag_index_map().size() == b.mesh_load().size());
            }
        }

        {
            bool all_empty = true;
            for (const Bin& b0: bin())
            {
                if (!b0.cell().empty())
                {
                    all_empty = false;
                    break;
                }
            }
            assert(!all_empty);
        }

        std::cout << "dur1: " << dur1 << " -- " << rank << std::endl;
        std::cout << "dur2: " << dur2 << " -- " << rank << std::endl;
        std::cout << "dur3: " << dur3 << " -- " << rank << std::endl;
        std::cout << "dur4: " << dur4 << " -- " << rank << std::endl;
        std::cout << "dur5: " << dur5 << " -- " << rank << std::endl;
    }

    void RegularMesh::register_resident_cell(const MeshCell& cell, const Tag& cell_tag, const Tag& mesh_tag, bool calcload, int rank)
    {
        bool overlap = aabb_.do_intersect(AABB(cell.poly()));
        if(!overlap)
        {
            return;
        }


        std::vector<BinRMTag> tag;
        get_bintag_adaptive(cell.poly().centroid(), tag, false); // centroid might be inbetween bins. so tags.size coudl be > 1.

        for (const BinRMTag& brmt: tag) 
        {
            //Bin& _bin = bin_p(tag.front()); // chose only one bin even if multiple tags exist.
            Bin& _bin = bin_p(brmt);

            assert(cell_tag.isvalid());
            assert(mesh_tag.isvalid());

            bool point_in_bin = _bin.is_resident(cell.poly().centroid());

            /*if (cell.tag()() == 18003 && mesh_tag() == 0)
            {
                std::cout << "bin min 0: " << _bin.aabb().min(0) << std::endl;
                std::cout << "bin min 1: " << _bin.aabb().min(1) << std::endl;
                std::cout << "bin min 2: " << _bin.aabb().min(2) << std::endl;

                std::cout << "bin max 0: " << _bin.aabb().max(0) << std::endl;
                std::cout << "bin max 1: " << _bin.aabb().max(1) << std::endl;
                std::cout << "bin max 2: " << _bin.aabb().max(2) << std::endl;

                std::cout << "cell cnt 0: " << cell.poly().centroid()(0) << std::endl;
                std::cout << "cell cnt 1: " << cell.poly().centroid()(1) << std::endl;
                std::cout << "cell cnt 2: " << cell.poly().centroid()(2) << std::endl;

                for (const auto& pp: cell.poly().vertices())
                {
                    std::cout << "cell vertex 0: " << pp.r(0) << std::endl;
                    std::cout << "cell vertex 1: " << pp.r(1) << std::endl;
                    std::cout << "cell vertex 2: " << pp.r(2) << std::endl;
                }

                assert(false);
            }*/

            if (point_in_bin)
            {
                _bin.add_bincell(cell_tag, mesh_tag, rank);
                assert(!_bin.cell().empty());
                if (calcload)
                {
                    _bin.increment_mesh_load(mesh_tag);
                }

                break;
            }
        }
    }

    void RegularMesh::register_cell(const MeshCell& cell, const Tag& cell_tag, const Tag& mesh_tag, bool calcload, int rank, bool adaptive, double& dur1, double& dur2, double& dur3, double& dur4, double& dur5)
    {
        //assert(!bintag_index_map_.left.empty());
        double start, stop;

        start = MPI_Wtime();
        bool overlap = aabb_.do_intersect(AABB(cell.poly()));
        stop = MPI_Wtime();
        dur1 += stop - start;
        if(!overlap)
        {
            return;
        }
        //assert(overlap);

        std::vector<BinRMTag> tag;
        tag.reserve(64);
        int dummy;
        for (const Bin& b: bin_)
        {
            assert(b.rm() == nullptr);
        }

        start = MPI_Wtime();
        get_bintag_adaptive(AABB(cell.poly()), tag, dummy);
        stop = MPI_Wtime();
        dur2 += stop - start;
        assert(!tag.empty());

        start = MPI_Wtime();
        std::vector<Tag> ut;
        ut.reserve(tag.size());
        for (const BinRMTag& _tag: tag)
            ut.push_back(_tag.rmtag());

        std::sort(ut.begin(), ut.end(), &Tag::sorter);
        ut.erase(std::unique(ut.begin(), ut.end()), ut.end());

        std::vector<RegularMeshIndex> minb(ut.size());
        std::vector<RegularMeshIndex> maxb(ut.size());

        assert(!ut.empty());
        stop = MPI_Wtime();
        dur3 += stop - start;
        
        start = MPI_Wtime();

        // determine corner bins.
        //RegularMeshIndex minb, maxb;
        for (int i=0; i<ut.size(); ++i)
            //for (const Tag& _ut: ut)
        {
            int mini = BIG_POS_NUM;
            int minj = BIG_POS_NUM;
            int mink = BIG_POS_NUM;
            int maxi = BIG_NEG_NUM;
            int maxj = BIG_NEG_NUM;
            int maxk = BIG_NEG_NUM;
            for (const BinRMTag& _tag: tag)
            {
                if (_tag.rmtag() == ut[i])
                {
                    if (adaptive)
                    {
                        const Bin& _bin = bin(_tag);
                        mini = std::min(mini, _bin.index().i());
                        minj = std::min(minj, _bin.index().j());
                        mink = std::min(mink, _bin.index().k());
                        maxi = std::max(maxi, _bin.index().i());
                        maxj = std::max(maxj, _bin.index().j());
                        maxk = std::max(maxk, _bin.index().k());
                    }
                    else
                    {
                        const Bin& _bin = bin(_tag.bintag());
                        mini = std::min(mini, _bin.index().i());
                        minj = std::min(minj, _bin.index().j());
                        mink = std::min(mink, _bin.index().k());
                        maxi = std::max(maxi, _bin.index().i());
                        maxj = std::max(maxj, _bin.index().j());
                        maxk = std::max(maxk, _bin.index().k());
                    }

                }
            }

            minb[i].set(mini, minj, mink);
            maxb[i].set(maxi, maxj, mink);
        }

        stop = MPI_Wtime();
        dur4 += stop - start;

        assert(!ut.empty());
        assert(!minb.empty());
        assert(!maxb.empty());

        start = MPI_Wtime();

        RegularMeshIndex ind;
        for (int u=0; u<ut.size(); ++u)
        //for (const BinRMTag& _tag: tag)
        {
            //for (const Tag& t: tag)
            for (int k=minb[u].k(); k<=maxb[u].k(); ++k)
            {
                for (int i=minb[u].i(); i<=maxb[u].i(); ++i)
                {
                    //assert(t.isvalid());
                    for (int j=minb[u].j(); j<=maxb[u].j(); ++j)
                    {
                        ind.set(i, j, k);
                        Bin& _bin = bin_p(ind.i(), ind.j(), ind.k());
                        //Bin& _bin = bin_p(_tag);

                        //assert(_bin.aabb().do_intersect(AABB(cell.polygon())));
                        if (_bin.aabb().do_intersect(AABB(cell.poly())))
                        {
                            //BinCell bc;
                            //bc.set_cell(cell_tag);
                            //bc.set_mesh(mesh_tag);
                            assert(cell_tag.isvalid());
                            assert(mesh_tag.isvalid());
                            //std::cout << "cell_tag = " << cell_tag() << std::endl;
                            //std::cout << "mesh_tag = " << mesh_tag() << std::endl;
                            //_bin.add_bincell(bc);
                            _bin.add_bincell(cell_tag, mesh_tag, rank);
                            assert(!_bin.cell().empty());

                            //if (!cell.is_ghost())
                            {
                                bool point_in_bin = _bin.is_resident(cell.poly().centroid());
                                //bool point_in_bin = point_inside_bin(ind, cell.polygon().centroid());
                                //bool point_in_bin = point_inside_bin(ind, cell.polytope()->centroid());
                                if (point_in_bin)
                                {
                                    if (calcload)
                                    {
                                        _bin.increment_mesh_load(mesh_tag);
                                        assert(_bin.mesh_tag_index_map().size() == _bin.mesh_load().size());
                                    }
                                }
                            }
                        }
                    }
                }
        }
        }

        stop = MPI_Wtime();
        dur5 += stop - start;
    }
}
