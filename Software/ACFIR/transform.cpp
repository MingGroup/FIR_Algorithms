
#include "transform.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace dsplib{
    namespace trans{
        int bitrev(int k, int len) {
            int kr = 0;
            for(int i = 0; i < len; i++) {
                kr |= ((k >> i) & 1) << (len - 1 - i);
            }
            return kr;
        }

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

        cvfp fft_forward(const cvfp& coef) {
            int n = 1, logn = 0;
            for(int i = (coef.size()-1); i > 0; i /= 2) {n *= 2; logn++;}
            cvfp eval(n);
            for(int i = 0; i < coef.size(); i++) eval[i] = coef[i];
            for(int i = coef.size(); i < n; i++) eval[i] = {0, 0};

            // Prepare rootList
            cvfp RootList_rev(n / 2);
            for(int i = 0; i < n / 2; i++) {
                RootList_rev[bitrev(i, logn - 1)] = {cos(2.0 * pi * i / n), -sin(2.0 * pi * i / n)};
            }
            for(int s = 0; s < logn; s++) {
                int m = 1 << s;
                int k = 1 << (logn - 1 - s); // the width of the bunch
                // std::cout << ">>>>> STAGE = " << s << std::endl;
                for(int i = 0; i < m; i++) {
                    // cnfp W = RootList_rev[m+i];
                    cnfp W = RootList_rev[i];
                    // std::cout << " - - ADDRESS = " << bitrev(i, logn-1) << std::endl;
                    // ui32 W = omegaList_rev[m+i];
                    for(int j = 0; j < k; j++) {
                        int j0 = 2 * k * i + j;
                        int j1 = 2 * k * i + j + k;
                        cnfp A = eval[j0];
                        cnfp B = eval[j1];
                        
                        cnfp T = B * W;
                        cnfp E = A + T;
                        cnfp O = A - T;
                        eval[j0] = E;
                        eval[j1] = O;
                    }
                }
            }
            return eval;
        }

        cvfp fft_backward(const cvfp& eval) {
            int n = 1, logn = 0;
            for(int i = (eval.size()-1); i > 0; i /= 2) {n *= 2; logn++;}
            cvfp coef(n);
            for(int i = 0; i < coef.size(); i++) coef[i] = eval[i];
            for(int i = coef.size(); i < n; i++) coef[i] = {0, 0};

            // Prepare rootList
            cvfp invRootList_rev(n / 2);
            for(int i = 0; i < n / 2; i++) {
                invRootList_rev[bitrev(i, logn - 1)] = {cos(2.0 * pi * i / n), sin(2.0 * pi * i / n)};
            }
            for(int s = 0; s < logn; s++) {
                int m = 1 << (logn - 1 - s);
                int k = 1 << s;
                for(int i = 0; i < m; i++) {
                    cnfp W = invRootList_rev[i];
                    // ui32 W = inv_omegaList_rev[m+i];
                    for(int j = 0; j < k; j++) {
                        int j0 = 2*k*i + j;
                        int j1 = 2*k*i + j + k;
                        cnfp A = coef[j0];
                        cnfp B = coef[j1];
                        cnfp E = A + B;
                        cnfp O = A - B;
                        O = O * W;
                        coef[j0] = E / 2.0;
                        coef[j1] = O / 2.0;
                    }
                }
            }
            return coef;
        }
    }
}