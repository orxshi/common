#include "remap.h"

namespace Common
{
    void Remapper::profstart(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->start(fun);
    }

    void Remapper::profstop(std::string fun)
    {
        if (profiler_ == nullptr) {
            return;
        }

        profiler_->stop(fun);
    }

    Remapper::Remapper(SpatialPartitionContainer* spc, const MPI_Comm& comm, std::shared_ptr<Profiler> profiler, bool verbose, bool remap_batch): spc_(spc), world_(comm, boost::mpi::comm_attach), profiler_(profiler), remap_batch_(remap_batch), verbose_(verbose)
    {
    }

    //void Remapper::remap(Loadmap& lm, bool greedy)
    void Remapper::remap(Loadmap& lm, bool mergebins)
    {
        if (verbose_) {
            std::cout << "inside remap function - " << world_.rank() << std::endl;
        }
        std::vector<MeshTransferInfo> mti;
        //profstart("second-remap-dist");
        if (verbose_) {
            std::cout << "distributing - " << world_.rank() << std::endl;
        }

        lm.distribute_mti_2(mti);
        //profstop("second-remap-dist");
        if (verbose_) {
            std::cout << "second distributed - " << world_.rank() << std::endl;
        }
        if (mti.empty())
        {
            std::cout << "mti is empty for rank " << world_.rank() << std::endl;
        }
        assert(!mti.empty());
        //profstart("spsender-cotr");
        //std::function<void(const SpatialPartition&)> f = [&](const SpatialPartition& s) {add_sp(s);};
        SpSender spsender(world_, profiler_);
        //profstop("spsender-cotr");
        //profstart("second-remap-send");
        //spsender.send(lm.mesh(), lm.rm(), mti, false);
        std::deque<SpatialPartition> incoming;

        if (verbose_) {
            std::cout << "spc spsender cotr - " << world_.rank() << std::endl;
        }

        profstart("spsender-transfer");
        if (remap_batch_) {
            spsender.transfer_batch(lm.mesh(), lm.rm(), mti, false, incoming);
        }
        else {
            spsender.transfer_per_dest(lm.mesh(), lm.rm(), mti, false, incoming);
        }
        for (const auto& ss: incoming)
        {
            assert(!ss.aabb().faces().empty());
        }
        profstop("spsender-transfer");

        if (verbose_) {
            std::cout << "transfer batch - " << world_.rank() << std::endl;
        }

        //profstart("spc-make-cell-tags");
        //make_cell_tags(incoming);
        //profstop("spc-make-cell-tags");

        std::vector<bool> cond(incoming.size(), false);

        if (verbose_) {
            std::cout << "remap just a test - " << world_.rank() << std::endl;
        }

        //if (!greedy)
        if (mergebins)
        {
            assert(false);
            if (verbose_) {
                std::cout << "reducing - " << world_.rank() << std::endl;
            }
            profstart("reduce");
            //spc_->reduce_sp(incoming);
            profstop("reduce");
            
            if (verbose_) {
                std::cout << "reduced - " << world_.rank() << std::endl;
            }
        }
        else
        {
            if (verbose_) {
                std::cout << "not reducing - " << world_.rank() << std::endl;
            }

            if (world_.rank() != 0)
            {
                assert(!incoming.empty());
            }

            profstart("spc-add-incoming");
            spc_->add_sp(incoming, cond);
            profstop("spc-add-incoming");

            if (verbose_) {
                std::cout << "added incoming sp - " << world_.rank() << std::endl;
            }

            profstart("spc-push-incoming");
            spc_->push_sp(incoming, cond);
            profstop("spc-push-incoming");

            if (verbose_) {
                std::cout << "pushed incoming sp - " << world_.rank() << std::endl;
            }

            if (world_.rank() != 0)
            {
                assert(!spc_->sp().empty());
            }
        }

        for (const auto& sp: spc_->sp())
        {
            assert(!sp.aabb().faces().empty());
        }

        if (verbose_) {
            std::cout << "removing dup cells - " << world_.rank() << std::endl;
        }

        //profstart("spc-merge-incoming");
        //merge_sp(incoming);
        //profstop("spc-merge-incoming");

        //if (!mergebins)
        //{
            //profstart("spc-remove-dup");
            //spc_->remove_dup_cells_and_points();
            //profstop("spc-remove-dup");
        //}

        if (verbose_) {
            std::cout << "remap test - " << world_.rank() << std::endl;
        }

        if (spc_ != nullptr) {
            if (verbose_) {
                std::cout << "asserting test - " << world_.rank() << std::endl;
            }
            assert(duplicate_exist(*spc_) == false);
            if (verbose_) {
                std::cout << "asserted test - " << world_.rank() << std::endl;
            }
        }

        //profstop("second-remap-send");
        if (verbose_) {
            std::cout << "second send - " << world_.rank() << std::endl;
        }
        //profstart("second-remap-remove");
        lm.remove_cells_from_rm();
        spc_->set_global_rm(lm.rm());
        spc_->make_mesh_aabbs();
        //profstop("second-remap-remove");
        assert(spc_->global_rm().tag()() == 0);

        //if (print_mesh_) {
            //print_all_meshes_in_partitions();
        //}
    }
}
