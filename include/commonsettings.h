#ifndef COMMON_SETTINGS_H
#define	COMMON_SETTINGS_H

#include <boost/program_options.hpp>
#include <boost/serialization/utility.hpp>
#include <fstream>
#include <iostream>
#include <boost/mpi.hpp>

namespace Common
{
    struct Settings
    {
        friend class boost::serialization::access;

        public:

        boost::program_options::variables_map vm_;

        virtual void config() = 0;

        private:
    };
}

#endif
