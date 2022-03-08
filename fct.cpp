/*
Radix-2 Out-of-Place DIT FFT Algorithm for 1D Real Input

TODO:
    - implement IDFT
    - replace complex number operations with macros
    - add SIMD build support
*/

#include <complex>
#include <math.h>
#include <iostream>

using std::cout;
using std::endl;
using std::conj;
using std::cos;
using std::sin;


typedef double real;
typedef std::complex<real> complex;


const complex i = complex(0.0, 1.0);
const double pi=3.141592653589793238462643383279502884197169399375105820974944;


/*
Calculates twiddle factors (complex roots of unity) for an N/2-point DFT
and the associated post-processing for a real-valued 1D input
*/
complex* calculate_fft_twiddles(int N) {
    complex* twiddles = new complex[N];
    double phase;
    for (int k=0; k<N; k++) {
        phase = -2.0 * pi * k / N;
        twiddles[k] = complex(cos(phase), sin(phase));
    }
    return twiddles;
}


/*
Combine the outputs of two DFTs
*/
void r2_butterfly(
    complex* twiddles,
    complex* output,
    int stride,
    int m
) {
    complex* output2 = output + m;
    complex t;
    do {
        t = *output2 * *twiddles;
        *output2 = *output - t;
        *output += t;

        twiddles += 2 * stride;
        output++;
        output2++;
    } while (--m);
}


/*
Radix-2 Cooley-Tukey FFT
*/
void fft_recursive(
    complex* twiddles,
    const complex* input,
    int n,
    complex* output,
    int stride
) {
    int m = n/2;

    if (m == 1) {
        output[0] = input[0];
        output[1] = input[stride];
    } else {
        fft_recursive(twiddles, input, m, output, 2*stride);
        fft_recursive(twiddles, input+stride, m, output+m, 2*stride);
    }

    r2_butterfly(twiddles, output, stride, m);
}


/*
Collapses real input of size N into complex sequence of size N/2,
calculates the N/2-point DFT, then extracts the N-point DFT of the
original N-point real input sequence
*/
void fft(
    complex* twiddles,
    const real* input,
    int N,
    complex* output
) {

    // collapse real input into N/2-point complex sequence
    complex* input_complex = new complex[N];
    cout << "N/2 COLLAPSED: ";
    for (int k=0; k<N/2; k++) {
        input_complex[k] = complex(input[2 * k], input[2 * k + 1]);
        cout << '|' << input_complex[k];
    }
    cout << "|\n";

    // perform N/2-point FFT on complex sequence
    fft_recursive(twiddles, input_complex, N/2, output, 1);
    cout << "N/2-POINT FFT OUTPUT: ";
    for (int k=0; k<N/2; k++) {
        cout << '|' << output[k];
    }
    cout << "|\n";

    // derive N-point DFT of input data from N/2-point DFT
    output[N/2] = output[0].real() - output[0].imag();
    output[0] = output[0].real() + output[0].imag();
    output[3*N/4] = output[N/4];
    output[N/4] = conj(output[N/4]);
    complex T, Tc, c3, c4, c5;
    for (int k=1; k<N/4; k++) {
        T = output[k];
        Tc = conj(output[N/2 - k]);

        c3 = T + Tc;
        c4 = T - Tc;
        c5 = i * twiddles[k] * c4;

        output[k] = 0.5 * (c3 - c5);
        output[N/2 - k] = conj(0.5 * (c3 + c5));
        output[N - k] = conj(output[k]);
        output[N/2 + k] = conj(output[N/2 - k]);

    }

}


complex* calculate_fct_twiddles(int N) {
    complex* twiddles = new complex[N];
    real phase;
    for (int k=0; k<N; k++) {
        phase = -pi * k / (2.0 * N);
        twiddles[k] = complex(cos(phase), sin(phase));
    }
    return twiddles;
}


void reorder(
    real* input,
    real* output,
    int N
) {
    int n = N >> 1;
    for (int k=0; k<n; k++) {
        output[k] = input[2 * k];
        output[n + k] = input[2 * k + 1];
    }
}


void fct(
    complex* fct_twiddles,
    real* input,
    real* output,
    complex* fft_output,
    complex* fft_twiddles,
    int N
) {
    reorder(input, output, N);
    cout << "REORDERED: ";
    for (int k=0; k<N; k++) {
        cout << '|' << output[k];
    }
    cout << "|\n";

    fft(fft_twiddles, output, N, fft_output);
    cout << "FFT OUTPUT: ";
    for (int k=0; k<N; k++) {
        cout << '|' << fft_output[k];
    }
    cout << "|\n";

    for (int k=0; k<N; k++) {
        output[k] = (fct_twiddles[k] * fft_output[k]).real();
    }
}


int main() {

    int N = 16;
    complex* fct_twiddles = calculate_fct_twiddles(N);
    real* input = new real[N];
    real* output = new real[N];
    complex* fft_output = new complex[N];
    complex* fft_twiddles = calculate_fft_twiddles(N);

    for (int k=0; k<N; k++) {
        input[k] = 1 + k%2 * -1;
    }

    cout << "FCT INPUT: ";
    for (int k=0; k<N; k++) {
        cout << '|' << input[k];
    }
    cout << "|\n";

    fct(fct_twiddles, input, output, fft_output, fft_twiddles, N);

    cout << "FCT OUTPUT: ";
    for (int k=0; k<N; k++) {
        cout << '|' << 2 * output[k];
    }
    cout << "|\n";

    return 0;
}
