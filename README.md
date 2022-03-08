# fast_dct

Fast Discrete Cosine Transform (FCT) using an FFT optimzed for real input sequences. Original paper by can be found [here](http://eelinux.ee.usm.maine.edu/courses/ele486/docs/makhoul.fastDCT.pdf) with the only change being that the algorithm implemented in `fct.cpp` does not multiply the output by a scale factor of two.

FFT used is a Radix-2 Cooley-Tukey algorithm optimized for real input sequences [here](https://github.com/hrichharms/cooley-tukey_fft).
