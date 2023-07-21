
#include "transform.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace dsplib{
    namespace trans{
        cvfp dft_forward(const cvfp& coef) {
            int n = 1;
            for(int i = (coef.size()-1); i > 0; i /= 2) n *= 2;
            cnfp Wn(cos(2.0 * pi / n), -sin(2.0 * pi / n));
            cvfp eval(n, {0, 0});
            for(int k = 0; k < n; k++) {
                cnfp W = {1, 0};
                cnfp Wd = pow(Wn, k);
                for(int j = 0; j < n; j++) {
                    cnfp value;
                    if(j >= coef.size()) {
                        value = {0, 0};
                    } else {
                        value = coef[j];
                    }
                    eval[k] += W * value;
                    W = W * Wd;
                }
            }
            return eval;
        }

        cvfp dft_backward(const cvfp& eval) {
            int n = 1;
            for(int i = (eval.size()-1); i > 0; i /= 2) n *= 2;
            cnfp Wn(cos(2.0 * pi / n), sin(2.0 * pi / n));
            cvfp coef(n, {0, 0});
            for(int k = 0; k < n; k++) {
                cnfp W = {1, 0};
                cnfp Wd = pow(Wn, k);
                for(int j = 0; j < n; j++) {
                    cnfp value;
                    if(j >= eval.size()) {
                        value = {0, 0};
                    } else {
                        value = eval[j];
                    }
                    coef[k] += W * value;
                    W = W * Wd;
                }
                coef[k] /= n;
            }
            return coef;
        }
    }
}