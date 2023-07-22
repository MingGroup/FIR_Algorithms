
#include "evaluation.h"

namespace dsplib{
    namespace eval{
        double L_2_error(const vfp& v0, const vfp& v1) {
            double err = 0;
            int len = (v0.size() > v1.size()) ? v0.size() : v1.size();
            for(int i = 0; i < len; i++) {
                double e0 = (i < v0.size()) ? v0[i] : 0.0;
                double e1 = (i < v1.size()) ? v1[i] : 0.0;
                err += (e0 - e1) * (e0 - e1);
            }
            err /= len;
            return err;
        }
    }
}