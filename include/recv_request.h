#ifndef COMMON_RECV_REQUEST_H
#define	COMMON_RECV_REQUEST_H

namespace Common
{
    struct RecvRequest
    {
        int source;
        BinRMTag sp_tag;
        bool received;
        boost::mpi::request request;
        
        RecvRequest(int s, const BinRMTag& t): source(s), sp_tag(t), received(false)
        {
        }

        RecvRequest(): received(false)
        {
        }
    };
}

#endif
