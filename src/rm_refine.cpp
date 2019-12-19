#include "regular_mesh.h"

namespace Common
{
    void RegularMesh::refine_adaptive(const std::deque<Mesh>& meshes, const BinRMTag& heaviest_bt, int rank, bool uniproc, int& new_rmtag, int& new_bt)
    {
        auto _rmtag = rmtag_address_map_.find(heaviest_bt.rmtag());
        assert(_rmtag != rmtag_address_map_.end());

        auto _rm = _rmtag->second;
        assert(_rm != nullptr);

        Bin& heavy_bin = _rm->bin_p(heaviest_bt.bintag());
        new_rmtag = get_new_rmtag();
        new_bt = get_new_bt();
        heavy_bin.init_rm(new_bt, new_rmtag, nstripe_);
        insert_to_rmtag_address_map(new_rmtag, heavy_bin.rm());
        if (uniproc || rank != 0) {
            heavy_bin.register_cells_to_rm(meshes);
        }
    }
    
    void RegularMesh::refine(const std::vector<Mesh>& meshes, int rank)
    {
        // Define new step length.
        vec3<int> new_nstripe(2 * nstripe(0), 2 * nstripe(1), 2 * nstripe(2));
        //unsigned int new_nstripe = 2. * nstripe();
        vec3<double> newh = h()/2.;
        // Create new bins.
        std::vector<Bin> newbin;
        newbin.resize(new_nstripe(0)*new_nstripe(1)*new_nstripe(2));

        // set indices of new bins.
        for (int depth=0; depth<nstripe(2); ++depth)
        {
            for (int row=0; row<nstripe(1); ++row)
            {
                for (int col=0; col<nstripe(0); ++col)
                {
                    // indices of lower left quadrant.
                    int ll_col = 2. * col;
                    int ll_row = 2. * row;
                    int ll_depth = 2. * depth;
                    newbin[ll_col + ll_row*new_nstripe(0) + ll_depth*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(ll_row, ll_col, ll_depth));

                    //
                    int col_ = ll_col + 1;
                    int row_ = ll_row;
                    int depth_ = ll_depth;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                    //
                    col_ = ll_col;
                    row_ = ll_row + 1;
                    depth_ = ll_depth;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                    //
                    col_ = ll_col + 1;
                    row_ = ll_row + 1;
                    depth_ = ll_depth;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                    ////

                    col_ = ll_col + 1;
                    row_ = ll_row;
                    depth_ = ll_depth + 1;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                    //
                    col_ = ll_col;
                    row_ = ll_row + 1;
                    depth_ = ll_depth + 1;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));

                    //
                    col_ = ll_col + 1;
                    row_ = ll_row + 1;
                    depth_ = ll_depth + 1;
                    newbin[col_ + row_*new_nstripe(0) + depth_*new_nstripe(0)*new_nstripe(1)].set_index(RegularMeshIndex(row_, col_, depth_));
                }
            }
        }

        // iterate through cells of unrefined bin and register them to new bins.
        // register cell to all overlapping bins in AABB.
        for (int depth=0; depth<nstripe(2); ++depth)
        {
            for (int row=0; row<nstripe(1); ++row)
            {
                for (int col=0; col<nstripe(0); ++col)
                {
                    RegularMeshIndex indorg(row, col, depth);
                    Bin& binorg = bin_p(indorg); // reference to original bin.
                    int ll_col = 2. * col;
                    int ll_row = 2. * row;
                    int ll_depth= 2. * depth;
                    int col_ = ll_col + 1;
                    int row_ = ll_row + 1;
                    int depth_ = ll_depth + 1;
                    RegularMeshIndex ind(ll_row, ll_col, ll_depth); // lower-left index of new bin.
                    // loop through new bins which are subset of unrefined bin.
                    for (int k=ll_depth; k<=depth_; ++k)
                    {
                        for (int i=ll_row; i<=row_; ++i)
                        {
                            for (int j=ll_col; j<=col_; ++j)
                            {
                                ind.set(i, j, k);
                                Bin& _bin = newbin[ind.j()+ind.i()*new_nstripe(0)+ind.k()*new_nstripe(0)*nstripe(1)]; // reference to new bin.
                                //_bin.resize_mesh_load(binorg.mesh_load().size());
                                // make a quad to represent new bin.
                                double d0 = aabb().min(0) + ind.j() * newh(0);
                                double d1 = aabb().min(1) + ind.i() * newh(1);
                                double d2 = aabb().min(2) + ind.k() * newh(2);
                                vec3<double> llc(d0, d1, d2);
                                vec3<double> ur = llc + newh;
                                _bin.set_aabb(AABB(llc, ur));
                                //std::vector<Point> pts = {Point(llc), Point(llc(0)+newh(0), llc(1)), Point(llc+newh), Point(llc(0), llc(1)+newh(1))};
                                //Polygon _quad(pts);
                                //AABB quad(_quad);
                                //Polytope quad(pts);
                                // loop through bincells of unrefined bin.
                                for (int t=0; t<binorg.cell().size(); ++t)
                                {
                                    Tag mt = binorg.cell(t).mesh();
                                    Tag ct = binorg.cell(t).cell();
                                    for (auto m=meshes.begin(); m!=meshes.end(); ++m)
                                    {
                                        if (mt == m->tag())
                                        {
                                            if (_bin.aabb().do_intersect(AABB(m->cell(ct).poly())))
                                                //if (quad.do_intersect(*m->cell(ct).polytope()))
                                            {
                                                //BinCell bc;
                                                //bc.set_cell(ct);
                                                //bc.set_mesh(mt);
                                                assert(ct.isvalid());
                                                assert(mt.isvalid());
                                                //_bin.add_bincell(bc);
                                                _bin.add_bincell(ct, mt, rank);

                                                //if (!m->cell(ct).is_ghost())
                                                {
                                                    //bool point_in_bin = quad.do_intersect(m->cell(ct).polygon().centroid());
                                                    bool point_in_bin = _bin.aabb().do_intersect(m->cell(ct).poly().centroid());
                                                    if (point_in_bin)
                                                    {
                                                        _bin.increment_mesh_load(mt);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                }
                }
            }
    }

        set_nstripe(new_nstripe);
        set_h(newh);
        bin_ = newbin;
        assert(bin_.size() == newbin.size());
    }

    /*void RegularMesh::refine(const std::vector<Mesh>& meshes)
      {
    // Define new step length.
    unsigned int new_nstripe = 2. * nstripe;
    vec3 newh = h/2.;
    // Create new bins.
    std::vector<Bin> newbin;
    newbin.resize(new_nstripe*new_nstripe, Bin(bin[0].mesh_load.size()));

    for (int row=0; row<nstripe; ++row)
    {
    for (int col=0; col<nstripe; ++col)
    {
    // indices of lower left quadrant.
    int ll_col = 2. * col;
    int ll_row = 2. * row;
    newbin[ll_col + ll_row*new_nstripe].row = ll_row;
    newbin[ll_col + ll_row*new_nstripe].col = ll_col;

    //
    int col_ = ll_col + 1;
    int row_ = ll_row;
    newbin[col_ + row_*new_nstripe].row = row_;
    newbin[col_ + row_*new_nstripe].col = col_;

    //
    col_ = ll_col;
    row_ = ll_row + 1;
    newbin[col_ + row_*new_nstripe].row = row_;
    newbin[col_ + row_*new_nstripe].col = col_;

    //
    col_ = ll_col + 1;
    row_ = ll_row + 1;
    newbin[col_ + row_*new_nstripe].row = row_;
    newbin[col_ + row_*new_nstripe].col = col_;
    }
    }

    // iterate through cells of unrefined bin and register them to new bins.
    // register cell to all overlapping bins in AABB.
    for (int row=0; row<nstripe; ++row)
    {
    for (int col=0; col<nstripe; ++col)
    {
    Index indorg(row, col, 0);
    Bin& binorg = (*this)(indorg);
    int ll_col = 2. * col;
    int ll_row = 2. * row;
    int col_ = ll_col + 1;
    int row_ = ll_row + 1;
    Index ind(ll_row, ll_col, 0);
    for (ind.i=ll_row; ind.i<=row_; ++ind.i)
    {
    for (ind.j=ll_col; ind.j<=col_; ++ind.j)
    {
    Bin& bin_ = newbin[ind.j+ind.i*new_nstripe];
    vec3 llc;
    llc[0] = aabb.min[0] + ind.j * newh[0];
    llc[1] = aabb.min[1] + ind.i * newh[1];
    std::vector<Point> pts = {Point(llc), Point(llc[0]+newh[0], llc[1]), Point(llc+newh), Point(llc[0], llc[1]+newh[1])};
    Polygon quad(pts);
    for (int t=0; t<binorg.cell.size(); ++t)
    {
    int mt = binorg.cell[t].mesh;
    int ct = binorg.cell[t].cell;
    if (quad.do_intersect(meshes[mt].cell[ct].polygon))
    {
    BinCell bc;
    bc.cell = ct;
    bc.mesh = mt;
    //bin_.cell_tag.push_back(ct);
    //bin_.mesh_tag.push_back(mt);
    // make point in bin intersection test to determine centroid owner.
    vec3 p = meshes[mt].cell[ct].polygon.centroid().r;
    bool point_in_bin = false;
    if (p[0] >= llc[0] && p[0] <= llc[0]+newh[0])
    {
        if (p[1] >= llc[1] && p[1] <= llc[1]+newh[1])
        {
            point_in_bin = true;
        }
    }
    //bin_.centroid_owner.push_back(point_in_bin);
    bc.resident = point_in_bin;
    bin_.cell.push_back(bc);
    if (point_in_bin)
    {
        ++bin_.mesh_load[mt];
    }
}
}
}
}
}
}

nstripe = new_nstripe;
h = newh;
bin = newbin;
}*/
}
