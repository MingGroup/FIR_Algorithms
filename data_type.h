
#ifndef DATA_TYPE_H
#define DATA_TYPE_H

#include <inttypes.h>
#include <vector>
#include <complex>
#include <cmath>

namespace dsplib{
    const double pi = acos(-1.0);
    typedef double nfp;
    typedef std::complex<double> cnfp;
    typedef std::vector<nfp>  vfp; // full precision vector
    typedef std::vector<cnfp> cvfp; // full precision complex vector

    cvfp vfp_2_cvfp(const vfp& vi, int len);
    std::vector<vfp> cvfp_2_vfp(const cvfp& vi, int len);
};

#endif