#ifndef COMMON_MOVING_CELL_H
#define	COMMON_MOVING_CELL_H

#include "commonmesh.h"

namespace Common
{
    struct MovingCellPre
    {
        std::vector<int> mesh_;
        std::vector<int> cell_;
        std::vector<unsigned int> dest_;
        unsigned int bin_;
        std::vector<const MeshCell*> address_;

        MovingCellPre(int cell, int mesh, unsigned int dest, unsigned int bin, const MeshCell* addr);
        void add(int cell, int mesh, unsigned int dest, const MeshCell* addr);
    };
    
    class MovingCell
    {
        std::vector<MovingCellPre> mcp_;
        typedef std::vector<std::vector<MeshCell>> cell_t;
        cell_t cell_;
        std::map<unsigned int, unsigned int> rank_to_index_;
        std::map<unsigned int, unsigned int> index_to_rank_;
        std::map<unsigned int, unsigned int> binindex_to_index_;
        std::map<unsigned int, unsigned int> index_to_binindex_;
        //std::vector<std::pair<Tag,Tag>> indices_;

        public:

        void resize_cell();
        const std::map<unsigned int, unsigned int>& index_to_binindex() const;
        void add_cell_real(const std::deque<Mesh>& mesh, std::string s);
        void add_cell_pre(unsigned int rank, const MeshCell& c, unsigned int binindex);
        void add_cell(unsigned int rank, const MeshCell& c, unsigned int binindex);
        unsigned int index_to_rank(unsigned int i);
        unsigned int index_to_binindex(unsigned int i);
        const cell_t& cell() const;
        const std::vector<MeshCell>& cell(unsigned int i) const;
        void clear();
        void reserve_for_mcp(size_t size);
    };
}

#endif
