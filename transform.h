
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "data_type.h"

namespace dsplib{
    namespace trans{
        cvfp dft_forward(const cvfp& coef);
        cvfp dft_backward(const cvfp& eval);
    }
}

#endif