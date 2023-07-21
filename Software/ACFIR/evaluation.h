
#ifndef EVALUATION_H
#define EVALUATION_H

#include "data_type.h"

namespace dsplib{
    namespace eval{
        /**
         * @brief Compute the L_2 error between 2 real vectors
         **/
        double L_2_error(const vfp& v0, const vfp& v1);
    }
}

#endif