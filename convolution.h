/*
    Convolution.h
    FIR is, to some extend, a vector convolution operator
*/

#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "data_type.h"

namespace dsplib{
    namespace conv{
        vfp fir_regular(const vfp& seq, const vfp& tap);
    }
}

#endif