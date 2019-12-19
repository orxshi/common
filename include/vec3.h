#ifndef COMMON_VEC2_H
#define	COMMON_VEC2_H

#include <mpi.h>
#include <boost/mpi.hpp>
#include <iostream>

namespace Common
{
    template<class T>
        class vec3
        {
            public:

            vec3(T x, T y, T z);
            vec3(): vec3(0., 0., 0.) {};

            const T& operator()(int i) const;
            void set(T x, T y, T z);
            void set_x(T v);
            void set_y(T v);
            void set_z(T v);
            void set(T v) = delete;
            vec3& operator=(const vec3& other);
            bool operator==(const vec3& other) const;
            bool operator!=(const vec3& other) const;
            vec3 operator-(const vec3& x) const;
            vec3 operator+(const vec3& x) const;
            vec3 operator/(T x) const;
            vec3 operator*(T x) const;
            T len();
            template<class Archive> void serialize(Archive & ar, const unsigned int version)
            {
                ar & data_;
            }

            private:

            T data_[3];

            friend class boost::serialization::access;
        };

    template<class T>
        vec3<T> normal(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c);
    template<class T>
        vec3<T> cross(const vec3<T>& a, const vec3<T>& b);
    template<class T>
        T dotp(const vec3<T>& a, const vec3<T>& b);
}

#include "vec3.hpp"

#endif
