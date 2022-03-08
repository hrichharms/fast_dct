# fast_dct

Fast Discrete Cosine Transform (FCT) using FFT. Original paper by can be found [here](http://eelinux.ee.usm.maine.edu/courses/ele486/docs/makhoul.fastDCT.pdf). Unlike the algorithm described in the aforementioned paper, the FCT implemented in `fct.cpp` does not multiply the outputs by a scale factor of two.

FFT used is a Radix-2 Cooley-Tukey algorithm optimized for real input sequences [here](https://github.com/hrichharms/cooley-tukey_fft).
