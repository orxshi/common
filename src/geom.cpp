#include "geom.h"

namespace Common
{
    int orientation (const vec3<double>& p, const vec3<double>& q, const vec3<double>& r)
    {
        // Find orientation of ordered triplet (p, q, r).
        // Return following values.
        // 0 --> p, q and r are colinear
        // 1 --> clockwise
        // 2 --> counterclockwise

        // See 10th slides from following link for derivation of the formula
        // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf

        double val = (q(1) - p(1)) * (r(0) - q(0)) - (q(0) - p(0)) * (r(1) - q(1));

        //if (val == 0.) return 0;  // colinear
        if (std::abs(val) < ZERO) return 0;  // colinear

        return (val > 0.) ? 1 : 2; // clock or counterclockwise
    }
}
