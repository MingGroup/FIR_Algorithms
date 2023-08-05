
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "transform.h"
#include "evaluation.h"
#include "convolution.h"

int main() {

    int n = 1024;
    // int n = 4096;
/*
    Part 0. FIR data Generation
*/
    // Generate input sequence
    dsplib::vfp xseq(n, 0);
    for(int i = 0; i < n; i++) {
        xseq[i] = (dsplib::nfp) (rand() % 2048);
    }

    // Assign tapping coefficients of FIR
    dsplib::vfp tap(n, 0);
    for(int i = 0; i < n; i++) {
        tap[i] = (dsplib::nfp) (rand() % 2048);
    }

/*
    Part 1. Direct FIR
    -> fir_o = xseq * tap
*/
    // FIR output
    auto fir_o = dsplib::conv::fir_regular(xseq, tap);
    std::cout << ">> FIR output: " << std::endl;
    for(int i = 0; i < fir_o.size(); i++) {
        std::cout << " | " << fir_o[i];
    }std::cout << std::endl;

/*
    Part 2. FIR by DFT
    -> dft_o = IDFT(DFT(xseq) .* DFT(tap))
*/
    // DFT(xseq)
    auto cxseq = dsplib::vfp_2_cvfp(xseq, 2*n);
    auto cxseq_t = dsplib::trans::dft_forward(cxseq);

    // DFT(tap)
    auto ctap = dsplib::vfp_2_cvfp(tap, 2*n);
    auto ctap_t = dsplib::trans::dft_forward(ctap);

    // point-wise multiplication between DFT(xseq) & DFT(tap)
    dsplib::cvfp mult_t(2*n);
    for(int i = 0; i < 2*n; i++) {
        mult_t[i] = cxseq_t[i] * ctap_t[i];
    }

    // IDFT(DFT(xseq) .* DFT(tap))
    auto cdft_o = dsplib::trans::dft_backward(mult_t);
    auto dft_o = dsplib::cvfp_2_vfp(cdft_o, 2*n);
    std::cout << ">> IDFT(DFT(xseq).*DFT(tap)): " << std::endl;
    for(int i = 0; i < dft_o[0].size(); i++) {
        std::cout << " | " << dft_o[0][i];
    }std::cout << std::endl;

/*
    Part 3. FIR by FFT
    -> fft_o = IFFT(DFT(xseq) .* FFT(tap))
*/
    // DFT(xseq)
    auto cxseq_ft = dsplib::trans::fft_forward(cxseq);

    // DFT(tap)
    auto ctap_ft = dsplib::trans::fft_forward(ctap);

    // point-wise multiplication between FFT(xseq) & FFT(tap)
    dsplib::cvfp mult_ft(2*n);
    for(int i = 0; i < 2*n; i++) {
        mult_ft[i] = cxseq_ft[i] * ctap_ft[i];
    }

    // IDFT(DFT(xseq) .* DFT(tap))
    auto cfft_o = dsplib::trans::fft_backward(mult_ft);
    auto fft_o = dsplib::cvfp_2_vfp(cfft_o, 2*n);
    std::cout << ">> IDFT(DFT(xseq).*DFT(tap)): " << std::endl;
    for(int i = 0; i < fft_o[0].size(); i++) {
        std::cout << " | " << fft_o[0][i];
    }std::cout << std::endl;

/*
    Part 4. Error Evaluation
*/
    // Error Evaluation
    double eval_L_2 = dsplib::eval::L_2_error(fir_o, dft_o[0]);
    std::cout << ">> The L_2 error between the output of DFT and FIR is : " << std::endl;
    std::cout << " | " << eval_L_2 << std::endl;
    eval_L_2 = dsplib::eval::L_2_error(fir_o, fft_o[0]);
    std::cout << ">> The L_2 error between the output of FFT and FIR is : " << std::endl;
    std::cout << " | " << eval_L_2 << std::endl;

    return 0;
}