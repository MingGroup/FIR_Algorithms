
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "data_type.h"

namespace dsplib{
    namespace trans{
        /**
         * @brief The Discrete Fourier Transform (DFT)
         * @param coef Input complex vector in Time Domain
         * @note The input will be automatically padding to length of 2^k 
         **/
        cvfp dft_forward(const cvfp& coef);
        /**
         * @brief The Inverse Discrete Fourier Transform (IDFT)
         * @param eval Input complex vector in Frequency Domain
         * @note The input will be automatically padding to length of 2^k 
         **/
        cvfp dft_backward(const cvfp& eval);
        /**
         * @brief The Fast Fourier Transform (FFT)
         * @note This function still needs further completion
         **/
        cvfp fft_forward(const cvfp& coef);
        /**
         * @brief The Inverse Fast Fourier Transform (IFFT)
         * @note This function still needs further completion
         **/
        cvfp fft_backward(const cvfp& eval);
    }
}

#endif