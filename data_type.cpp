
#include <vector>
#include <complex>
#include <cmath>
#include "data_type.h"

namespace dsplib{

    cvfp vfp_2_cvfp(const vfp& vi, int len) {
        cvfp vo(len, {0, 0});
        for(int i = 0; i < len; i++) {
            if(i < vi.size()) {
                vo[i].real(vi[i]);
                vo[i].imag(0);
            }
        }
        return vo;
    }

    std::vector<vfp> cvfp_2_vfp(const cvfp& vi, int len) {
        std::vector<vfp> vo(2);
        for(int i = 0; i < len; i++) {
            if(i < vi.size()) {
                vo[0].push_back(vi[i].real());
                vo[1].push_back(vi[i].imag());
            } else {
                vo[0].push_back(0);
                vo[1].push_back(0);
            }
        }
        return vo;
    }
}
