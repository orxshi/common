#include "adt.h"

namespace Common
{
    void ADTPoint::rotate(double angle, int axis, const vec3<double>& rot_point)
    {
        RotationMatrix rm;

        {
            vec3<double> z;
            z.set(dim_[0]-rot_point(0), dim_[2]-rot_point(1), dim_[4]-rot_point(2));
            auto newz = rm.rotate(angle, axis, z);

            dim_[0] = newz(0)+rot_point(0);
            dim_[2] = newz(1)+rot_point(1);
            dim_[4] = newz(2)+rot_point(2);
        }
        {
            vec3<double> z;
            z.set(dim_[1]-rot_point(0), dim_[3]-rot_point(1), dim_[5]-rot_point(2));
            auto newz = rm.rotate(angle, axis, z);

            dim_[1] = newz(0)+rot_point(0);
            dim_[3] = newz(1)+rot_point(1);
            dim_[5] = newz(2)+rot_point(2);
        }
    }

    void ADTPoint::move(const vec3<double>& v)
    {
        assert(dim_.size() == 4);

        dim_[0] += v(0);
        dim_[1] += v(0);
        dim_[2] += v(1);
        dim_[3] += v(1);
    }

    void ADTPoint::set_idx(int i)
    {
        idx_ = i;
    }


    ADTPoint::ADTPoint(const vec3<double>& p, int idx): idx_(idx)
    {
        dim_.resize(ADT_VAR);
        for (int i=0; i<ADT_DIM; ++i)
        {
            dim_[i*2]   = BIG_POS_NUM;
            dim_[i*2+1] = BIG_NEG_NUM;
        }

        for (int i=0; i<ADT_DIM; ++i)
        {
            dim_[i*2]   = std::min(dim_[i*2]  , p(i));
            dim_[i*2+1] = std::max(dim_[i*2+1], p(i));
        }
    }

    const std::vector<double>& ADTPoint::dim() const
    {
        return dim_;
    }

    double ADTPoint::dim(int i) const
    {
        assert(i >= 0);
        assert(i < dim_.size());

        return dim_[i];
    }

    int ADTPoint::idx() const
    {
        return idx_;
    }

    bool ADTPoint::overlap(const ADTPoint& other, bool verbose) const
    {
        //if (other.idx() == 537)
        //verbose = true;
        //else
        //verbose = false;

        for (int i=0; i<ADT_DIM; ++i)
        {
            if (verbose)
            {
                std::cout << "idx_ = " << idx_ << std::endl;

                std::cout << "dim_[0] = " << dim_[0] << std::endl;
                std::cout << "dim_[1] = " << dim_[1] << std::endl;
                std::cout << "dim_[2] = " << dim_[2] << std::endl;
                std::cout << "dim_[3] = " << dim_[3] << std::endl;

                std::cout << "other.dim(0) = " << other.dim(0) << std::endl;
                std::cout << "other.dim(1) = " << other.dim(1) << std::endl;
                std::cout << "other.dim(2) = " << other.dim(2) << std::endl;
                std::cout << "other.dim(3) = " << other.dim(3) << std::endl;
            }
            double dis1 = dim_[i*2] - other.dim(i*2+1);
            double dis2 = dim_[i*2+1] - other.dim(i*2);
            //if (dim_[i*2] > other.dim(i*2+1)) return false;
            if (dis1 > ZERO) return false;
            //if (dim_[i*2+1] < other.dim(i*2)) return false;
            if (dis2 < ZERO) return false;
        }

        return true;
    }
}
