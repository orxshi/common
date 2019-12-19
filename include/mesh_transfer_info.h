#ifndef COMMON_MESH_TRANSFER_INFO_H
#define	COMMON_MESH_TRANSFER_INFO_H

#include "boost/serialization/list.hpp"
#include "commonmesh.h"

namespace Common
{
    class MeshTransferInfo
    {
        public:

        struct MTIMesh
        {
            Tag tag_;
            //std::list<Tag> cell_;
            std::vector<Tag> cell_;
            friend class boost::serialization::access;

            MTIMesh() = default;
            MTIMesh(const Tag& tag, const Tag& ctag);
            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & tag_;
                ar & cell_;
            }
        };

        struct Dest
        {
            int rank_;
            BinRMTag sp_tag_;
            std::vector<MTIMesh> mesh_;
            std::deque<Mesh> make_meshes(const std::deque<Mesh>& m, const AABB& aabb, bool check_residency) const;
            std::deque<Mesh> make_meshes_without_pts(const std::deque<Mesh>& m) const;
            std::vector<MeshCell> make_cells(const std::deque<Mesh>& m, const AABB& aabb, bool check_residency) const;
            friend class boost::serialization::access;

            Dest() = default;
            Dest(int rank, const Tag& mtag, const Tag& ctag, const BinRMTag& sp_tag);
            void add(const Tag& mtag, const Tag& ctag);
            int ncell() const;
            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & rank_;
                ar & mesh_;
                ar & sp_tag_;
            }
        };

        private:

        int source_;
        std::vector<Dest> dest_;
        std::vector<int> rank_;
        friend class boost::serialization::access;

        public:

        MeshTransferInfo() = default;
        MeshTransferInfo(int source);
        MeshTransferInfo(int source, int dest, const Tag& mtag, const Tag& ctag, const BinRMTag& sp_tag);

        void add(int dest, const Tag& mtag, const Tag& ctag, const BinRMTag& sp_tag);
        int source() const;
        const std::vector<Dest>& dest() const;
        bool bin_exist(const BinRMTag& sp_tag) const;
        bool dest_exist(int rank) const;
        const std::vector<int> rank() const;
        template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & source_;
            ar & dest_;
        }
    };

    int pairing_function(int a, int b);
}

#endif
