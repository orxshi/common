#include "rotationmatrix.h"

namespace Common
{
    void RotationMatrix::set_rx(double yaw)
    {
        rx[0] = 1.;
        rx[1] = 0.;
        rx[2] = 0.;

        rx[3] = 0.;
        rx[4] = std::cos(yaw);
        rx[5] = -std::sin(yaw);

        rx[6] = 0.;
        rx[7] = std::sin(yaw);
        rx[8] = std::cos(yaw);
    }

    void RotationMatrix::set_ry(double pitch)
    {
        ry[0] = std::cos(pitch);
        ry[1] = 0.;
        ry[2] = std::sin(pitch);

        ry[3] = 0.;
        ry[4] = 1.;
        ry[5] = 0.;

        ry[6] = -std::sin(pitch);
        ry[7] = 0.;
        ry[8] = std::cos(pitch);
    }

    void RotationMatrix::set_rz(double roll)
    {
        rz[0] = std::cos(roll);
        rz[1] = -std::sin(roll);
        rz[2] = 0.;

        rz[3] = std::sin(roll);
        rz[4] = std::cos(roll);
        rz[5] = 0.;

        rz[6] = 0.;
        rz[7] = 0.;
        rz[8] = 1.;
    }

    vec3<double> RotationMatrix::rotate(double angle, int axis, const vec3<double>& v)
    {
        if (axis == 0) {
            set_rx(angle);
            return mul(rx, v);
        }
        else if (axis == 1) {
            set_ry(angle);
            return mul(ry, v);
        }
        else if (axis == 2) {
            set_rz(angle);
            return mul(rz, v);
        }
    }

    vec3<double> RotationMatrix::mul(const trimat& m, const vec3<double>& v)
    {
        vec3<double> c;

        c.set_x(m[0] * v(0) + m[1] * v(1) + m[2] * v(2));
        c.set_y(m[3] * v(0) + m[4] * v(1) + m[5] * v(2));
        c.set_z(m[6] * v(0) + m[7] * v(1) + m[8] * v(2));

        return c;
    }
}
