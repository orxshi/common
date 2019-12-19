#ifndef COMMON_SP_CELL_MOVER_H
#define COMMON_SP_CELL_MOVER_H

#include "commonmesh.h"
#include "regular_mesh.h"
#include "movingcell.h"
#include "sp.h"

namespace Common
{
    class ArrivalInfo
    {
        public:
            void add_size(size_t s);
            void add_disp(int begin, int end);
            size_t size() const;
            const std::vector<std::pair<int,int>>& disp() const;

            ArrivalInfo();

        private:
            size_t narrival_;
            std::vector<std::pair<int,int>> disp_;
    };

    class SpCellMover
    {
        public:
            void complete_receive(int* send_proc_2, std::vector<MeshCell>& arrivals_);
            void try_to_receive(int* send_proc_2, std::vector<MeshCell>& arrivals_);
            void inform_receivers(MPI_Win& send_proc_win, int* send_proc_2, MPI_Win& send_size_win, int* send_size, std::vector<int>& receive_count, std::vector<int> receive_size);
            void reserve_moving_elements(const std::map<BinRMTag, int>& bintag_proc_map, const RegularMesh& global_rm_, size_t& ut_size, size_t& move_size, int& nlevel, size_t& preresi, size_t& postresi);
            //void reserve_moving_elements(const std::vector<Outline>& outline);
            void push_moving_elements();
            void send_moving_cells(int* send_proc_2, std::vector<MeshCell>& arrivals_);
            void send_moving_cells(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win, int* send_proc_2, MPI_Win& send_size_win, int* send_size);
            void send_moving_dirichlet(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win_dirichlet_, int* send_proc_dirichlet_, std::vector<int>& receive_count);
            void send_moving_walls(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win_wall_, int* send_proc_wall_, std::vector<int>& receive_count);
            void send_moving_empty(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win_empty_, int* send_proc_empty_, std::vector<int>& receive_count);
            void send_moving_farfield(std::vector<boost::mpi::request>& request, MPI_Win& send_proc_win_farfield_, int* send_proc_farfield_, std::vector<int>& receive_count);
            void recv_moving_cells_self();
            void recv_moving_walls_self();
            void recv_moving_farfield_self();
            void recv_moving_empty_self();
            void recv_moving_dirichlet_self();
            void recv_moving_cells(int source);
            void reserve_arrival_cells();
            void recv_moving_cells(const std::vector<MeshCell>& arrivals);
            void recv_wall_cells(int source, int tag);
            void recv_dirichlet_cells(int source, int tag);
            void recv_empty_cells(int source, int tag);
            void recv_farfield_cells(int source, int tag);
            const MovingCell& movcel() const;
            std::vector<MeshCell>& arrival_cell();
            std::vector<MeshCell>& arrival_wall_cell();
            std::vector<MeshCell>& arrival_dirichlet_cell();
            std::vector<MeshCell>& arrival_empty_cell();
            std::vector<MeshCell>& arrival_farfield_cell();
            void clear_arrival_cell();
            void clear_arrival_wall_cell();
            void clear_arrival_dirichlet_cell();
            void clear_arrival_empty_cell();
            void clear_arrival_farfield_cell();
            const Tag& tag() const;
            void add_arrival_info(size_t size, int disp_begin, int disp_end);
            void prepare_moving_elements(std::vector<SpCellMover>& sp_cell_mover);

            SpCellMover(SpatialPartition* sp, const MPI_Comm& comm, bool mergebins);

        private:
            std::vector<bool> rece;
            std::vector<MeshCell> arrival_cell_; 
            std::vector<MeshCell> arrival_wall_cell_; 
            std::vector<MeshCell> arrival_dirichlet_cell_; 
            std::vector<MeshCell> arrival_empty_cell_; 
            std::vector<MeshCell> arrival_farfield_cell_; 
            MovingCell movcel_;
            MovingCell wall_cell_;
            MovingCell dirichlet_cell_;
            MovingCell empty_cell_;
            MovingCell farfield_cell_;
            SpatialPartition* sp_;
            boost::mpi::communicator world_;
            Tag tag_;
            ArrivalInfo arrival_info_;
            bool mergebins_;

            bool master() const;
    };
}

#endif
