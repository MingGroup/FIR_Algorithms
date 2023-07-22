
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "transform.h"
#include "evaluation.h"
#include "convolution.h"

int main() {

    int n = 8;

    // Generate input sequence
    dsplib::vfp xseq(n, 0);
    for(int i = 0; i < n; i++) {
        xseq[i] = i;
    }

    // Assign tapping coefficients of FIR
    dsplib::vfp tap(n, 0);
    for(int i = 0; i < n; i++) {
        tap[i] = i;
    }

    // FIR output
    auto fir_o = dsplib::conv::fir_regular(xseq, tap);
    std::cout << ">> FIR output: " << std::endl;
    for(int i = 0; i < fir_o.size(); i++) {
        std::cout << " | " << fir_o[i];
    }std::cout << std::endl;

    // DFT(xseq)
    auto cxseq = dsplib::vfp_2_cvfp(xseq, 2*n);
    auto cxseq_t = dsplib::trans::dft_forward(cxseq);

    // DFT(tap)
    auto ctap = dsplib::vfp_2_cvfp(tap, 2*n);
    auto ctap_t = dsplib::trans::dft_forward(ctap);

    // point-wise multiplication between DFT(xseq) & DFT(tap)
    dsplib::cvfp mult(2*n);
    for(int i = 0; i < 2*n; i++) {
        mult[i] = cxseq_t[i] * ctap_t[i];
    }

    // IDFT(DFT(xseq) .* DFT(tap))
    auto dft_o_t = dsplib::trans::dft_backward(mult);
    auto dft_o = dsplib::cvfp_2_vfp(dft_o_t, 2*n);
    std::cout << ">> IDFT(DFT(xseq).*DFT(tap)): " << std::endl;
    for(int i = 0; i < dft_o[0].size(); i++) {
        std::cout << " | " << dft_o[0][i];
    }std::cout << std::endl;
    
    // Error Evaluation
    double eval_L_2 = dsplib::eval::L_2_error(fir_o, dft_o[0]);
    std::cout << ">> The L_2 error between the output of DFT and FIR is : " << std::endl;
    std::cout << " | " << eval_L_2 << std::endl;

    return 0;
}