#include "load_calculator.h"

namespace Common
{
    bool LoadCalculator::master() const
    {
        if (world_.rank() == 0)
        {
            return true;
        }

        return false;
    }

    LoadCalculator::LoadCalculator(LoadEstimType type, std::string area_rep, bool merge_bins, const MPI_Comm& comm): type_(type), world_(comm, boost::mpi::comm_attach), area_rep_(area_rep), merge_bins_(merge_bins)
    {
    }

    /*double LoadCalculator::area_load_convex(const Bin& bin, const std::deque<Mesh>& mesh) const
    {
        assert(!mesh.empty());

        double load = 0;
        if (bin.cell().empty())
        {
            return load;
        }

        if (bin.mesh_tag_index_map().left.size() == 1)
        {
            //load = mesh_load[ii];
            return load;
        }

        std::deque<std::vector<Point>> pts = bin.mesh_pts(mesh);
        assert(pts.size() == mesh.size());
        
        assert(!pts.empty());
        std::deque<Polygon_2> pols;
        for (const auto& pp: pts)
        {
            if (pp.empty())
            {
                pols.push_back(Polygon_2());
            }
            else
            {
                std::vector<CGAL_Point> ptscgal;
                for (const Point& p: pp)
                {
                    ptscgal.push_back(CGAL_Point(p.r(0), p.r(1)));
                }
                std::vector<CGAL_Point> outpts(ptscgal.size());
                auto loc = CGAL::convex_hull_2(ptscgal.begin(), ptscgal.end(), outpts.data());
                pols.push_back(Polygon_2(outpts.begin(), outpts.begin() + (loc - outpts.data())));
            }
        }

        const auto& mesh_load = bin.mesh_load();

        for (auto z=bin.mesh_tag_index_map().left.begin(); z!=bin.mesh_tag_index_map().left.end(); ++z)
        {
            int mt = z->first;
            int ii = z->second;

            auto mm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mt;});
            assert(mm != mesh.end());

            int di = std::distance(mesh.begin(), mm);

            if (std::next(z) == bin.mesh_tag_index_map().left.end())
            {
                break;
            }
            
            for (auto zz=std::next(z); zz!=bin.mesh_tag_index_map().left.end(); ++zz)
            {
                int mtt = zz->first;
                int jj = zz->second;

                auto mmm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mtt;});
                assert(mmm != mesh.end());

                int dj = std::distance(mesh.begin(), mmm);

                double oa = area_of_overlap(pols[di], pols[dj]);

                assert(oa >= 0.);

                //if (world_.rank() == 1)
                //{
                    //std::cout << "xs: " << pols[di].size() << std::endl;
                    //std::cout << "xx: " << CGAL::to_double(pols[di].vertex(0).x()) << std::endl;
                    //std::cout << "xy: " << CGAL::to_double(pols[di].vertex(0).y()) << std::endl;
                //}
                //
                double polareai = std::abs(CGAL::to_double(pols[di].area()));
                double polareaj = std::abs(CGAL::to_double(pols[dj].area()));

                double overratioi = oa / polareai;
                double overratioj = oa / polareaj;

                double nonoverratioi = 1. - overratioi;
                double nonoverratioj = 1. - overratioj;

                double alpha = 2.0;
                double beta = 1.0;

                double ncommon_i = mesh_load[ii] * (alpha * overratioi + beta * nonoverratioi);
                double ncommon_j = mesh_load[jj] * (alpha * overratioj + beta * nonoverratioj);

                //double ncommon_i = (oa * mesh_load[ii] / std::abs(CGAL::to_double(pols[di].area())));
                //double ncommon_j = (oa * mesh_load[jj] / std::abs(CGAL::to_double(pols[dj].area())));

                load += (ncommon_i + ncommon_j);
            }
        }

        return load;
    }*/

    /*double LoadCalculator::area_load_alpha(const Bin& bin, const std::deque<Mesh>& mesh) const
    {
        assert(!mesh.empty());

        double load = 0;
        if (bin.cell().empty())
        {
            return load;
        }

        std::deque<std::vector<Point>> pts = bin.mesh_pts(mesh);

        assert(!pts.empty());
        std::deque<Polygon_2> pols;
        for (const auto& pp: pts)
        {
            if (pp.empty())
            {
                pols.push_back(Polygon_2());
            }
            else
            {
                std::vector<CGAL_Point> ptscgal;
                for (const Point& p: pp)
                {
                    ptscgal.push_back(CGAL_Point(p.r(0), p.r(1)));
                }
                Alpha_shape_2 A(ptscgal.begin(), ptscgal.end(), FT(0.01), Alpha_shape_2::REGULARIZED);

                Alpha_shape_2::Alpha_shape_vertices_iterator vit = A.alpha_shape_vertices_begin(),
                    vend = A.alpha_shape_vertices_end();

                Alpha_shape_2::Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
                    end = A.alpha_shape_edges_end();

                //Triangulation_2::Edge_circulator ec = incident_edges(A.point(*vit));
                Triangulation_2::Vertex_circulator vc = A.incident_vertices(*vit);
                Triangulation_2::Vertex_circulator vccopy = A.incident_vertices(*vit);
                --vccopy;

                assert(vit != vend);
                std::vector<CGAL_Point> outpts;
                outpts.push_back(A.point(*vit));
                do
                {
                    //std::cout << "vc: " << vc->point() << " " << world_.rank() << std::endl;
                    outpts.push_back(vc->point());
                    //std::cout << "out: " << outpts.back() << " " << world_.rank() << std::endl;
                }
                while(++vc != vccopy);
                assert(!outpts.empty());

                //for(; it!=end; ++it)
                //{
                    //outpts.push_back(A.segment(*it).point(0));
                //}
                    
                for(const auto& pi: outpts)
                {
                    std::cout << "out - " << world_.rank() << " " << pi << std::endl;
                }

                //std::cout << "outsizeeeee: " << outpts.size() << std::endl;

                //std::vector<CGAL_Point> ppts = {CGAL_Point(-3.44, -10), CGAL_Point(-3.4068, -10.02), CGAL_Point(-3.36672, -10.02), CGAL_Point(-3.36, -10), CGAL_Point(0.5, 0.5)};
                //std::vector<CGAL_Point> ppts = outpts;
                //std::vector<CGAL_Point> ppts = {CGAL_Point(-1.35135, -31.1812), CGAL_Point(-1.35135, -31.0811), CGAL_Point(-1.45145, -30.981), CGAL_Point(-1.65165, -30.6807), CGAL_Point(-1.45145, -31.1812)};
                //Polygon_2 pp(ppts.begin(), ppts.end());

                pols.push_back(Polygon_2(outpts.begin(), outpts.end()));
                //std::cout << "cheking empty - " << world_.rank() << std::endl;
                //assert(!pols.back().is_empty());
                //std::cout << "cheking simple - " << world_.rank() << std::endl;
                //std::cout << "the first pol - " << world_.rank() << " " << pp << std::endl;
                //std::cout << "the reall pol - " << world_.rank() << " " << pols.back() << std::endl;
                //assert(pp.is_simple());
                //std::cout << "now huh - " << world_.rank() << std::endl;
                //assert(Polygon_2(outpts.begin(), outpts.end()).is_simple());
                //std::cout << "now lolo - " << world_.rank() << std::endl;
                //pols.back().is_simple();
                //std::cout << "now checl - " << world_.rank() << std::endl;
                std::cout << "checkiing - " << world_.rank() << std::endl;
                assert(pols.back().is_simple());
                std::cout << "check finished - " << world_.rank() << std::endl;

                //pols.push_back(PS::simplify(Polygon_2(outpts.begin(), outpts.end()), Cost(), StopBelow(50)));
            }

            //for (int ii=0; ii<pols.size(); ++ii)
            //{
                //const auto& pol = pols[ii];

            //}
        }

        auto mesh_load = bin.mesh_load();

        for (auto z=bin.mesh_tag_index_map().left.begin(); z!=bin.mesh_tag_index_map().left.end(); ++z)
        {
            int mt = z->first;
            int ii = z->second;

            auto mm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mt;});
            assert(mm != mesh.end());

            int di = std::distance(mesh.begin(), mm);

            if (std::next(z) == bin.mesh_tag_index_map().left.end()) break;
            
            for (auto zz=std::next(z); zz!=bin.mesh_tag_index_map().left.end(); ++zz)
            {
                int mtt = zz->first;
                int jj = zz->second;

                auto mmm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mtt;});
                assert(mmm != mesh.end());

                int dj = std::distance(mesh.begin(), mmm);

                double oa = area_of_overlap(pols[di], pols[dj]);

                assert(oa >= 0.);

                double ncommon_i = (oa * mesh_load[ii] / std::abs(CGAL::to_double(pols[di].area())));
                double ncommon_j = (oa * mesh_load[jj] / std::abs(CGAL::to_double(pols[dj].area())));

                load += (ncommon_i + ncommon_j);
            }
        }

        //std::cout << "load: " << load << std::endl;
        return load;
    }*/

    double LoadCalculator::area_load(const Bin& bin, const std::deque<Mesh>& mesh) const
    {
        assert(!mesh.empty());

        double load = 0;
        if (bin.cell().empty())
        {
            return load;
        }

        std::deque<AABB> aabbs = bin.mesh_aabb(mesh);
        assert(!aabbs.empty());
        for (const auto& aa: aabbs)
        {
            assert(!aa.faces().empty());
        }
        auto mesh_load = bin.mesh_load();

        for (auto z=bin.mesh_tag_index_map().left.begin(); z!=bin.mesh_tag_index_map().left.end(); ++z)
        {
            int mt = z->first;
            int ii = z->second;

            auto mm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mt;});
            assert(mm != mesh.end());

            int di = std::distance(mesh.begin(), mm);

            if (std::next(z) == bin.mesh_tag_index_map().left.end()) break;
            
            for (auto zz=std::next(z); zz!=bin.mesh_tag_index_map().left.end(); ++zz)
            {
                int mtt = zz->first;
                int jj = zz->second;

                auto mmm = std::find_if(mesh.begin(), mesh.end(), [&](const Mesh& m){return m.tag()() == mtt;});
                assert(mmm != mesh.end());

                int dj = std::distance(mesh.begin(), mmm);

                assert(dj < aabbs.size());
                double oa = aabbs[di].volume_of_overlap(aabbs[dj]);

                assert(oa >= 0.);
                //assert((oa / aabbs[di].area()) <= 1.);
                //assert((oa / aabbs[di].area()) >= 0.);
                /*if ((oa / aabbs[dj].area()) > 1.)
                {
                    std::cout << "oa: " << oa << std::endl;
                    std::cout << "area: " << aabbs[dj].area() << std::endl;
                    std::cout << "di: " << di << std::endl;
                    std::cout << "dj: " << dj << std::endl;
                    std::cout << "mm tag: " << mm->tag()() << std::endl;
                    std::cout << "mmm tag: " << mmm->tag()() << std::endl;
                    std::cout << "min0i: " << aabbs[di].min(0) << std::endl;
                    std::cout << "min1i: " << aabbs[di].min(1) << std::endl;
                    std::cout << "max0i: " << aabbs[di].max(0) << std::endl;
                    std::cout << "max1i: " << aabbs[di].max(1) << std::endl;
                    std::cout << "min0j: " << aabbs[dj].min(0) << std::endl;
                    std::cout << "min1j: " << aabbs[dj].min(1) << std::endl;
                    std::cout << "max0j: " << aabbs[dj].max(0) << std::endl;
                    std::cout << "max1j: " << aabbs[dj].max(1) << std::endl;
                }*/
                //assert((oa / aabbs[dj].area()) <= 1.);
                //assert((oa / aabbs[dj].area()) >= 0.);

                assert(!aabbs[di].faces().empty());
                assert(!aabbs[dj].faces().empty());
                double ncommon_i = (oa * mesh_load[ii] / aabbs[di].volume());
                double ncommon_j = (oa * mesh_load[jj] / aabbs[dj].volume());
                //double ncommon_i = oa * mesh_load[ii];
                //double ncommon_j = oa * mesh_load[jj];

                load += (ncommon_i + ncommon_j);
            }
        }

        return load;
    }

    double LoadCalculator::minmesh_load(const Bin& bin) const
    {
        if (bin.mesh_load().size() <= 1) {
            return 0.;
        }
        else {
            double load = 0.;
            auto ml = bin.mesh_load();
            for (auto it=ml.begin(); it!=ml.end(); ++it)
            {
                if (std::next(it) == ml.end()) {
                    break;
                }
                for (auto it2=std::next(it); it2!=ml.end(); ++it2)
                {
                    load += std::min(*it, *it2);
                }
            }
            return load;
        }
    }

    double LoadCalculator::solver_load(const Bin& bin) const
    {
        return std::accumulate(bin.mesh_load().begin(), bin.mesh_load().end(), 0.);
    }

    double LoadCalculator::hybrid_load(const Bin& bin) const
    {
        if (bin.mesh_load().size() <= 1) {
            return 0.;
        }
        else {
            return std::accumulate(bin.mesh_load().begin(), bin.mesh_load().end(), 0.);
        }
    }

    std::vector<double> mesh_load(const SpatialPartition& sp)
    {
        std::vector<double> mesh_load;
        for (const Mesh& m: sp.mesh())
        {
            double load = 0.;
            for (const MeshCell& mc: m.cell())
            {
                if (sp.is_resident(mc.poly().centroid()))
                {
                    ++load;
                }
            }

            mesh_load.push_back(load);
        }
        assert(!mesh_load.empty());
        return mesh_load;
    }

    /*std::vector<double> mesh_load(const SpatialPartition& sp, const Outline& outline)
    {
        std::vector<double> mesh_load;
        for (const Mesh& m: sp.mesh())
        {
            double load = 0.;
            for (const MeshCell& mc: m.cell())
            {
                //if (sp.is_resident(mc.polygon().centroid()))
                if (outline.do_contain(mc.polygon().centroid(), false))
                {
                    ++load;
                }
            }

            mesh_load.push_back(load);
        }
        assert(!mesh_load.empty());
        return mesh_load;
    }*/

    double LoadCalculator::load(const SpatialPartition& sp) const
    {
        if (type_ == LoadEstimType::hybrid)
        {
            return hybrid_load(sp);
        }
        else if (type_ == LoadEstimType::solver)
        {
            return solver_load(sp);
        }
        else if (type_ == LoadEstimType::minmesh)
        {
            assert(false);
            //return minmesh_load(sp, outline);
        }
        else if (type_ == LoadEstimType::area)
        {
            if (area_rep_ == "concave") {
                assert(false);
                //return area_load_alpha(sp);
            }
            else if (area_rep_ == "convex") {
                assert(false);
                //return area_load_convex(sp);
            }
            else if (area_rep_ == "aabb") {
                return area_load(sp);
            }
            else
            {
                assert(false);
            }
        }
        else
        {
            assert(false);
            return -1;
        }
    }

    /*double LoadCalculator::load(const SpatialPartition& sp, const Outline& outline) const
    {
        if (type_ == LoadEstimType::hybrid)
        {
            return hybrid_load(sp, outline);
        }
        else if (type_ == LoadEstimType::solver)
        {
            return solver_load(sp, outline);
        }
        else if (type_ == LoadEstimType::minmesh)
        {
            assert(false);
            //return minmesh_load(sp, outline);
        }
        else if (type_ == LoadEstimType::area)
        {
            if (settings_->arearep_ == "concave") {
                return area_load_alpha(sp, outline);
            }
            if (settings_->arearep_ == "convex") {
                return area_load_convex(sp, outline);
            }
            if (settings_->arearep_ == "aabb") {
                return area_load(sp, outline);
            }
        }
        else
        {
            assert(false);
            return -1;
        }
    }*/

    double LoadCalculator::load_r(const Bin& bin, const std::deque<Mesh>& mesh) const
    {
        if (type_ == LoadEstimType::hybrid)
        {
            return solver_load(bin);
        }
        else if (type_ == LoadEstimType::solver)
        {
            return solver_load(bin);
        }
        else if (type_ == LoadEstimType::minmesh)
        {
            assert(false);
            //return minmesh_load(bin);
        }
        else if (type_ == LoadEstimType::area)
        {
            if (area_rep_ == "concave") {
                assert(false);
                //return area_load_alpha(bin, mesh);
            }
            if (area_rep_ == "convex") {
                assert(false);
                //return area_load_convex(bin, mesh);
            }
            if (area_rep_ == "aabb") {
                return area_load(bin, mesh);
            }
        }
        else
        {
            assert(false);
            return -1;
        }
    }

    double LoadCalculator::load(const Bin& bin, const std::deque<Mesh>& mesh) const
    {
        if (type_ == LoadEstimType::hybrid)
        {
            return hybrid_load(bin);
        }
        else if (type_ == LoadEstimType::solver)
        {
            return solver_load(bin);
        }
        else if (type_ == LoadEstimType::minmesh)
        {
            assert(false);
            //return minmesh_load(bin);
        }
        else if (type_ == LoadEstimType::area)
        {
            if (area_rep_ == "concave") {
                assert(false);
                //return area_load_alpha(bin, mesh);
            }
            if (area_rep_ == "convex") {
                assert(false);
                //return area_load_convex(bin, mesh);
            }
            if (area_rep_ == "aabb") {
                return area_load(bin, mesh);
            }
        }
        else
        {
            assert(false);
            return -1;
        }
    }

    /*double LoadCalculator::minmesh_load(const SpatialPartition& sp, const Outline& outline) const
    {
        auto ml = mesh_load(sp, outline);

        if (ml.size() <= 1) {
            return 0.;
        }
        else {
            double load = 0.;
            for (auto it=ml.begin(); it!=ml.end(); ++it)
            {
                if (std::next(it) == ml.end()) {
                    break;
                }
                for (auto it2=std::next(it); it2!=ml.end(); ++it2)
                {
                    load += std::min(*it, *it2);
                }
            }
            return load;
        }
    }*/

    double LoadCalculator::hybrid_load(const SpatialPartition& sp) const
    {
        auto ml = mesh_load(sp);

        if (ml.size() <= 1) {
            return 0.;
        }
        else {
            return std::accumulate(ml.begin(), ml.end(), 0.);
        }
    }

    /*double LoadCalculator::hybrid_load(const SpatialPartition& sp, const Outline& outline) const
    {
        auto ml = mesh_load(sp, outline);

        if (ml.size() <= 1) {
            return 0.;
        }
        else {
            return std::accumulate(ml.begin(), ml.end(), 0.);
        }
    }*/
    
    /*double LoadCalculator::solver_load(const SpatialPartition& sp, const Outline& outline) const
    {
        auto ml = mesh_load(sp, outline);
        return std::accumulate(ml.begin(), ml.end(), 0.);
    }*/

    double LoadCalculator::solver_load(const SpatialPartition& sp) const
    {
        auto ml = mesh_load(sp);
        return std::accumulate(ml.begin(), ml.end(), 0.);
    }

    /*double LoadCalculator::area_load_convex(const SpatialPartition& sp) const
    {
        auto ml = mesh_load(sp);

        double load = 0;

        if (sp.mesh().empty()) {
            return load;
        }

        std::deque<std::vector<Point>> pts = sp.mesh_pts();
        if (!pts.empty());
        std::deque<Polygon_2> pols;
        for (const auto& pp: pts)
        {
            if (pp.empty())
            {
                pols.push_back(Polygon_2());
            }
            else
            {
                std::vector<CGAL_Point> ptscgal;
                for (const Point& p: pp)
                {
                    ptscgal.push_back(CGAL_Point(p.r(0), p.r(1)));
                }
                std::vector<CGAL_Point> outpts(ptscgal.size());
                auto loc = CGAL::convex_hull_2(ptscgal.begin(), ptscgal.end(), outpts.data());
                pols.push_back(Polygon_2(outpts.begin(), outpts.begin() + (loc - outpts.data())));
            }
        }

        std::deque<AABB> aabbs = sp.mesh_aabb();

        for (auto imesh = sp.mesh().begin(); imesh != sp.mesh().end(); ++imesh)
        {
            int ii = std::distance(sp.mesh().begin(), imesh);

            if (std::next(imesh) == sp.mesh().end()) break;

            for (auto jmesh = std::next(imesh); jmesh != sp.mesh().end(); ++jmesh)
            {
                int jj = std::distance(sp.mesh().begin(), jmesh);

                double oa = area_of_overlap(pols[ii], pols[jj]);

                assert(oa >= 0.);

                load += (oa * ml[ii] / CGAL::to_double(pols[ii].area())) + (oa * ml[jj] / std::abs(CGAL::to_double(pols[jj].area())));
            }
        }

        return load;
    }*/

    /*double LoadCalculator::area_load_convex(const SpatialPartition& sp, const Outline& outline) const
    {
        auto ml = mesh_load(sp, outline);

        double load = 0;

        if (sp.mesh().empty()) {
            return load;
        }

        std::deque<std::vector<Point>> pts = sp.mesh_pts(outline);
        if (!pts.empty());
        std::deque<Polygon_2> pols;
        for (const auto& pp: pts)
        {
            if (pp.empty())
            {
                pols.push_back(Polygon_2());
            }
            else
            {
                std::vector<CGAL_Point> ptscgal;
                for (const Point& p: pp)
                {
                    ptscgal.push_back(CGAL_Point(p.r(0), p.r(1)));
                }
                std::vector<CGAL_Point> outpts(ptscgal.size());
                auto loc = CGAL::convex_hull_2(ptscgal.begin(), ptscgal.end(), outpts.data());
                pols.push_back(Polygon_2(outpts.begin(), outpts.begin() + (loc - outpts.data())));
            }
        }

        std::deque<AABB> aabbs = sp.mesh_aabb(outline);

        for (auto imesh = sp.mesh().begin(); imesh != sp.mesh().end(); ++imesh)
        {
            int ii = std::distance(sp.mesh().begin(), imesh);

            if (std::next(imesh) == sp.mesh().end()) break;

            for (auto jmesh = std::next(imesh); jmesh != sp.mesh().end(); ++jmesh)
            {
                int jj = std::distance(sp.mesh().begin(), jmesh);

                double oa = area_of_overlap(pols[ii], pols[jj]);

                assert(oa >= 0.);

                load += (oa * ml[ii] / CGAL::to_double(pols[ii].area())) + (oa * ml[jj] / std::abs(CGAL::to_double(pols[jj].area())));
            }
        }

        return load;
    }*/

    /*double LoadCalculator::area_load_alpha(const SpatialPartition& sp) const
    {
        auto ml = mesh_load(sp);

        double load = 0;

        if (sp.mesh().empty()) {
            return load;
        }

        std::deque<std::vector<Point>> pts = sp.mesh_pts();
        if (!pts.empty());
        std::deque<Polygon_2> pols;
        for (const auto& pp: pts)
        {
            if (pp.empty())
            {
                pols.push_back(Polygon_2());
            }
            else
            {
                std::vector<CGAL_Point> ptscgal;
                for (const Point& p: pp)
                {
                    ptscgal.push_back(CGAL_Point(p.r(0), p.r(1)));
                }
                Alpha_shape_2 A(ptscgal.begin(), ptscgal.end(), FT(0.01), Alpha_shape_2::REGULARIZED);

                Alpha_shape_2::Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
                    end = A.alpha_shape_edges_end();

                std::vector<CGAL_Point> outpts;
                for(; it!=end; ++it)
                {
                    outpts.push_back(A.segment(*it).point(0));
                }

                //Alpha_shape_2::Alpha_shape_vertices_iterator it = A.alpha_shape_vertices_begin(),
                    //end = A.alpha_shape_vertices_end();

                //std::vector<CGAL_Point> outpts;
                //for(; it!=end; ++it)
                //{
                    //outpts.push_back(A.point(*it));
                //}

                pols.push_back(Polygon_2(outpts.begin(), outpts.end()));
            }
        }

        std::deque<AABB> aabbs = sp.mesh_aabb();

        for (auto imesh = sp.mesh().begin(); imesh != sp.mesh().end(); ++imesh)
        {
            int ii = std::distance(sp.mesh().begin(), imesh);

            if (std::next(imesh) == sp.mesh().end()) break;

            for (auto jmesh = std::next(imesh); jmesh != sp.mesh().end(); ++jmesh)
            {
                int jj = std::distance(sp.mesh().begin(), jmesh);

                double oa = area_of_overlap(pols[ii], pols[jj]);

                assert(oa >= 0.);

                load += (oa * ml[ii] / CGAL::to_double(pols[ii].area())) + (oa * ml[jj] / std::abs(CGAL::to_double(pols[jj].area())));
            }
        }

        return load;
    }*/

    /*double LoadCalculator::area_load_alpha(const SpatialPartition& sp, const Outline& outline) const
    {
        auto ml = mesh_load(sp, outline);

        double load = 0;

        if (sp.mesh().empty()) {
            return load;
        }

        std::deque<std::vector<Point>> pts = sp.mesh_pts(outline);
        if (!pts.empty());
        std::deque<Polygon_2> pols;
        for (const auto& pp: pts)
        {
            if (pp.empty())
            {
                pols.push_back(Polygon_2());
            }
            else
            {
                std::vector<CGAL_Point> ptscgal;
                for (const Point& p: pp)
                {
                    ptscgal.push_back(CGAL_Point(p.r(0), p.r(1)));
                }
                Alpha_shape_2 A(ptscgal.begin(), ptscgal.end(), FT(0.01), Alpha_shape_2::REGULARIZED);

                Alpha_shape_2::Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
                    end = A.alpha_shape_edges_end();

                std::vector<CGAL_Point> outpts;
                for(; it!=end; ++it)
                {
                    outpts.push_back(A.segment(*it).point(0));
                }

                //Alpha_shape_2::Alpha_shape_vertices_iterator it = A.alpha_shape_vertices_begin(),
                    //end = A.alpha_shape_vertices_end();

                //std::vector<CGAL_Point> outpts;
                //for(; it!=end; ++it)
                //{
                    //outpts.push_back(A.point(*it));
                //}

                pols.push_back(Polygon_2(outpts.begin(), outpts.end()));
            }
        }

        std::deque<AABB> aabbs = sp.mesh_aabb(outline);

        for (auto imesh = sp.mesh().begin(); imesh != sp.mesh().end(); ++imesh)
        {
            int ii = std::distance(sp.mesh().begin(), imesh);

            if (std::next(imesh) == sp.mesh().end()) break;

            for (auto jmesh = std::next(imesh); jmesh != sp.mesh().end(); ++jmesh)
            {
                int jj = std::distance(sp.mesh().begin(), jmesh);

                double oa = pols[ii].volume_of_overlap(pols[jj]);

                assert(oa >= 0.);

                load += (oa * ml[ii] / CGAL::to_double(pols[ii].area())) + (oa * ml[jj] / std::abs(CGAL::to_double(pols[jj].area())));
            }
        }

        return load;
    }*/

    double LoadCalculator::area_load(const SpatialPartition& sp) const
    {
        auto ml = mesh_load(sp);

        double load = 0;

        if (sp.mesh().empty()) {
            return load;
        }

        std::deque<AABB> aabbs = sp.mesh_aabb();
        assert(!aabbs.empty());

        for (auto imesh = sp.mesh().begin(); imesh != sp.mesh().end(); ++imesh)
        {
            int ii = std::distance(sp.mesh().begin(), imesh);

            if (std::next(imesh) == sp.mesh().end()) break;

            for (auto jmesh = std::next(imesh); jmesh != sp.mesh().end(); ++jmesh)
            {
                int jj = std::distance(sp.mesh().begin(), jmesh);

                assert(jj < aabbs.size());
                double oa = aabbs[ii].volume_of_overlap(aabbs[jj]);

                assert(oa >= 0.);
                //assert((oa / aabbs[ii].area()) <= 1.);
                //assert((oa / aabbs[ii].area()) >= 0.);
                //assert((oa / aabbs[jj].area()) <= 1.);
                //assert((oa / aabbs[jj].area()) >= 0.);

                assert(!aabbs[ii].faces().empty());
                assert(!aabbs[jj].faces().empty());
                load += (oa * ml[ii] / aabbs[ii].volume()) + (oa * ml[jj] / aabbs[jj].volume());
                //load += oa * (ml[ii] + ml[jj]);
            }
        }

        return load;
    }

    /*double LoadCalculator::area_load(const SpatialPartition& sp, const Outline& outline) const
    {
        auto ml = mesh_load(sp, outline);

        double load = 0;

        if (sp.mesh().empty()) {
            return load;
        }

        std::deque<AABB> aabbs = sp.mesh_aabb(outline);

        for (auto imesh = sp.mesh().begin(); imesh != sp.mesh().end(); ++imesh)
        {
            int ii = std::distance(sp.mesh().begin(), imesh);

            if (std::next(imesh) == sp.mesh().end()) break;

            for (auto jmesh = std::next(imesh); jmesh != sp.mesh().end(); ++jmesh)
            {
                int jj = std::distance(sp.mesh().begin(), jmesh);

                double oa = aabbs[ii].area_of_overlap(aabbs[jj]);

                assert(oa >= 0.);
                //assert((oa / aabbs[ii].area()) <= 1.);
                //assert((oa / aabbs[ii].area()) >= 0.);
                //assert((oa / aabbs[jj].area()) <= 1.);
                //assert((oa / aabbs[jj].area()) >= 0.);

                load += (oa * ml[ii] / aabbs[ii].area()) + (oa * ml[jj] / aabbs[jj].area());
                //load += oa * (ml[ii] + ml[jj]);
            }
        }

        return load;
    }*/

    bool LoadCalculator::load(const SpatialPartitionContainer& spc, double refine_tol, int iteration) const
    {
        bool remake_loadmap = false;
        //load_ = 0;
        double tload = 0;
        std::vector<double> _loads;

        if (!master())
        {
            if (merge_bins_)
            {
                assert(false);
                //auto outline = std::find_if(spc.outline().begin(), spc.outline().end(), [&](const Outline& ol){return ol.tag()() == world_.rank();});
                //for (const SpatialPartition& s: spc.sp())
                //{
                    //tload += load(s, *outline);
                //}
            }
            else
            {
                for (const SpatialPartition& s: spc.sp())
                {
                    tload += load(s);
                }
            }
        }

        gather(world_, tload, _loads, MASTER);

        double dev;
        double dev_ave;
        if (master())
        {
            auto result = std::minmax_element(_loads.begin()+1, _loads.end());
            double minload = *result.first;
            double maxload = *result.second;
            //for (auto ii=_loads.begin()+1; ii!=_loads.end(); ++ii) {
                //std::cout << "lcoad: " << *ii << std::endl;
            //}

            double average = std::accumulate(_loads.begin()+1, _loads.end(), 0.) / (_loads.size() - 1);
            dev_ave = std::max(maxload-average, average-minload);
            dev_ave *= 100. / average;

            dev = (maxload - minload) * 100. / maxload;
            //dev_ = dev;
            if (dev > refine_tol)
            {
                remake_loadmap = true;
            }

            std::fstream out;
            out.open("dev-vs-iter.csv", std::fstream::out | std::fstream::app);
            if (iteration == 0)
                out << "iter,dev\n";
            out << iteration << "," << dev << "\n";
            out.close();

            out.open("devave-vs-iter.csv", std::fstream::out | std::fstream::app);
            if (iteration == 0)
                out << "iter,dev\n";
            out << iteration << "," << dev_ave << "\n";
            out.close();

            out.open("load_iter.csv", std::fstream::out | std::fstream::app);
            if (iteration == 0)
            {
                out << "iter,";
                for (int i=1; i<world_.size(); ++i)
                {
                    out << "proc" << i;
                    if (i != world_.size() - 1)
                        out << ",";
                }
                out << "\n";
            }
            out << iteration << ",";
            for (int i=1; i<world_.size(); ++i)
            {
                out << _loads[i];
                if (i != world_.size() - 1)
                    out << ",";
            }
            out << "\n";
            out.close();
        }

        broadcast(world_, remake_loadmap, MASTER);

        return remake_loadmap;
    }

    void LoadCalculator::gather_load(const RegularMesh& rm, std::vector<double>& load_local, std::vector<double>& load_local_r, int& j, const std::deque<Mesh>& mesh)
    {
        for (const Bin& bin: rm.bin())
        {
            if (bin.rm() == nullptr)
            {
                //std::cout << "rank - " << world_.rank() << " bin nullptr" << std::endl;
                load_local[j] = load(bin, mesh);
                load_local_r[j] = load_r(bin, mesh);
                ++j;
            }
            else
            {
                //std::cout << "rank - " << world_.rank() << " bin not nullptr" << std::endl;
                gather_load(*bin.rm(), load_local, load_local_r, j, mesh);
            }
        }
    }

    /*void gather_load(const RegularMesh& rm, std::vector<double>& load_global, const std::deque<Mesh>& mesh, LoadEstimType load_estim_type)
    {
        size_t rm_size = size();
        load_global.resize(rm_size, 0);

        int j = 0;
        gather_load(rm, load_global, j, mesh, load_estim_type);
    }*/

    void LoadCalculator::gather_load(const RegularMesh& rm, std::vector<double>& load_global, std::vector<double>& load_global_r, const std::deque<Mesh>& mesh)
    {
        size_t rm_size = rm.size();
        assert(rm_size != 0);
        std::vector<double> load_local(rm_size, 0.);
        std::vector<double> load_local_r(rm_size, 0.);
        load_global.resize(rm_size, 0.);
        load_global_r.resize(rm_size, 0.);
        if (!master())
        {
            int j = 0;
            gather_load(rm, load_local, load_local_r, j, mesh);
        }
        //if (world_.rank() != 0)
        //{
            //auto localmax = *std::minmax_element(load_local.begin(), load_local.end()).second;
            //std::cout << "rank - " << world_.rank() << " loadlocal size " << load_local.size() << " localmax " << localmax << std::endl;
            ////std::cout << "rank - " << world_.rank() << " rm.bin size " << rm.bin().size() << std::endl;
            //assert(localmax != 0);
        //}
        //for (const auto& aa: load_local)
        //{
            //std::cout << "rank - " << world_.rank() << " load " << aa <<  std::endl;
        //}
        //

        reduce(world_, load_local.data(), load_local.size(), load_global.data(), std::plus<double>(), 0.);
        reduce(world_, load_local_r.data(), load_local_r.size(), load_global_r.data(), std::plus<double>(), 0.);

        /*if (load_estim_type_ == LoadEstimType::hybrid)
        {
            std::vector<int> dep(load_global.size(), -1);

            if (master())
            {
                std::vector<BinRMTag> tag;
                rm_->get_bin_tags(tag);

                for (int j=0; j<load_global.size(); ++j)
                {
                    if (load_global[j] != 0.) {
                        continue;
                    }

                    Tag rmtag = tag[j].rmtag();
                    Tag bintag = tag[j].bintag();

                    auto myrm = rm_->rmtag_address_map().find(rmtag);
                    if (myrm == rm_->rmtag_address_map().end())
                    {
                        std::cout << "rmtag: " << rmtag() << std::endl;
                        std::cout << "rmtag size: " << rm_->rmtag_address_map().size() << std::endl;
                        std::cout << "rm front: " << rm_->rmtag_address_map().begin()->first() << std::endl;
                        std::cout << "root rm tag: " << rm_->tag()() << std::endl;
                    }
                    assert(myrm != rm_->rmtag_address_map().end());

                    for (const Bin& b: myrm->second->bin())
                    {
                        if (b.rm() != nullptr)
                        {
                            // need to examine relative positions of sub-bins.
                        }
                        const Tag& t = b.tag();

                        if (t == bintag) {
                            continue;
                        }

                        auto it = std::find_if(tag.begin(), tag.end(), [&](const BinRMTag& tt){return tt == BinRMTag(t, rmtag);});
                        if (it == tag.end())
                        {
                            std::cout << "rmtag: " << rmtag() << std::endl;
                            std::cout << "t: " << t() << std::endl;
                            std::cout << "tag size: " << tag.size() << std::endl;
                            for (int z=0; z<tag.size(); ++z) 
                            {
                                std::cout << "tagg: " << tag[z].rmtag()() << " " << tag[z].bintag()() << std::endl;
                            }
                        }
                        assert(it != tag.end());
                        int dis = std::distance(tag.begin(), it);

                        if (load_global[dis] == 0.) {
                            continue;
                        }

                        assert(load_global_r[j] != 0);
                        load_global[dis] += load_global_r[j];
                        dep[j] = dis;
                    }
                }

                for (int i=0; i<dep.size(); ++i)
                {
                    if (dep[i] != -1)
                    {
                        int part = -1;
                        for (int uu=0; uu<aug_bin_tag_.size(); ++uu)
                        {
                            const auto& u = aug_bin_tag_[uu];
                            auto it = std::find_if(u.begin(), u.end(), [&](const BinRMTag& brmt){return brmt == tag[dep[i]];});
                            if (it != u.end()) {
                                part = uu;
                                break;
                            }
                        }

                        aug_bin_tag_[part].push_back(tag[i]);
                    }
                }
            }
        }*/
    }

    std::vector<BinRMTag> LoadCalculator::sorted_bin_tags(const RegularMesh& rm, std::vector<double>& load_global, std::vector<double>& load_global_r, const std::deque<Mesh>& mesh)
    {
        /*bool uniproc = false;

        if (comm.size() == 1) {
            uniproc = true;
        }*/

        /*if (uniproc) {
            gather_load(load_global, mesh, load_estim_type);
        }
        else {*/
            //std::cout << "rank - " << world_.rank() << " gathering load" << std::endl;
            gather_load(rm, load_global, load_global_r, mesh);
        //}

        std::vector<BinRMTag> tag;
        std::vector<BinRMTag> tag2;

        if (master())
        {
            rm.get_bin_tags(tag);
            assert(!tag.empty());
            tag2 = tag;
            assert(tag.size() == load_global.size());
            std::vector<int> index(tag.size());
            std::iota(index.begin(), index.end(), 0);
            std::sort(std::begin(index), std::end(index),
                    [&](int t1, int t2)
                    {
                    return load_global[t1] > load_global[t2];
                    });

            for (int i=0; i<tag.size(); ++i)
            {
                tag2[i] = tag[index[i]];
            }
        }

        if (master()) {
            assert(!tag2.empty());
        }
        return tag2;
    }
}
