#ifndef COMMON_ROTATIONMATRIX_H
#define	COMMON_ROTATIONMATRIX_H

#include <vec3.h>

namespace Common
{
    class RotationMatrix
    {
        public: 

        vec3<double> rotate(double angle, int axis, const vec3<double>& v);

        private:

        using trimat = std::array<double, 9>;
        trimat rx, ry, rz;

        void set_rx(double yaw);
        void set_ry(double pitch);
        void set_rz(double roll);
        vec3<double> mul(const trimat& m, const vec3<double>& v);
    };
}

#endif
