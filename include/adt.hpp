#ifndef COMMON_ADT_HPP
#define	COMMON_ADT_HPP

namespace Common
{
    template<typename Rai> ADTPoint::ADTPoint(Rai begin, Rai end, int idx): idx_(idx)
    {
        dim_.resize(ADT_VAR);
        for (int i=0; i<ADT_DIM; ++i)
        {
            dim_[i*2]   = BIG_POS_NUM;
            dim_[i*2+1] = BIG_NEG_NUM;
        }

        for (auto it=begin; it!=end; ++it)
        {
            for (int i=0; i<ADT_DIM; ++i)
            {
                dim_[i*2]   = std::min(dim_[i*2]  , it->r(i));
                dim_[i*2+1] = std::max(dim_[i*2+1], it->r(i));
            }
        }
    }
}

#endif
