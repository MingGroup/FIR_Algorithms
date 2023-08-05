

#include "convolution.h"

namespace dsplib{
    namespace conv{
        vfp fir_regular(const vfp& seq, const vfp& tap) {
            int slen = seq.size();
            int tlen = tap.size();
            int olen = slen + tlen - 1;
            vfp vo(olen);
            for(int i = 0; i < olen; i++) {
                vo[i] = 0;
                for(int t = 0; t <= i; t++) {
                    if((t < tlen) && ((i - t) < slen))
                        vo[i] += tap[t] * seq[i - t];
                }
            }
            return vo;
        };
    }
}