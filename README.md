# Signal-Separation


## Background Statistics Folder
- `whitening.m`: script which generates plot of whitening transformation on Gaussian random variables
- `joint_distributions.m`: script which generates plots of joint distributions of Gaussian and non-Gaussian random variables before and after orthogonal transformation
- `cross_corr_ex.m`: script that generates figures showing an example of finding time lag with cross correlation


## ICA Folder
This folder contains scripts and functions used to test independent component analysis, specifically fastica.
### Functions
- `preprocess.m`: Function that whitens and centers data
    - returns X: original data
    - X_mean: mean of original data
- `FP_kurt.m`: Fixed point iteration kurtosis ICA algorithm
    - returns S: estimated independent components
    - W: mixing matrix
    - w1: values of 1st row of mixing matrix W for each iteration of gradient descent
    - w2: values of 2nd row of mixing matrix W for each iteration of gradient descent
- `FP_negen.m`: Fixed point iteration negnetropy ICA algorithm
    - returns S: estimated independent components
    - W: mixing matrix
    - w1: values of 1st row of mixing matrix W for each iteration of gradient descent
    - w2: values of 2nd row of mixing matrix W for each iteration of gradient descent
- `calc_fft.m`: calcualtes the FFT of a given input signal (used for reordering signals in analysis and evaluation)
    - returns freq: frequencies from FFT
    - S_fft: spectrum (amplitudes) from FFT

### Comparison of Algorithms
- `sine_nonoise.m`
    - runs 1000 iterations for example of sine wave with 100 and 40 hz synthetic mixture without noise
    - compares fixed point iterations for kurtosis and negentropy
- `sine_noise.m`
    - runs 1000 iterations for example of sine wave with 100 and 40 hz synthetic mixture with noise
    - compares fixed point iterations for kurtosis and negentropy
- `k_wave_sim_no_timelag`
    - runs 1000 iterations for example from kwave (60 hz and 25 hz) 
    - does not account for time lag
- `k_wave_sim_timelag`
    - runs 1000 iterations for example from kwave (60 hz and 23 hz) 
    - accounts for time lag by direct computation from distance between source/sensor and speed of sound in medium defined
- `real_data_simulation`
    - runs 1000 iterations for example with real test data loaded in
- `optimization_plots.m`
    - Plots mixed component signals
    - Plots independent components after separation (normalized)
    - Plots optimization landscape of kurtosis and/or negentropy as abjective functions  
    - Plots FFT of independent components after separation
    - Plots FFT of original signals before mixing
    - Plots FFT of mixed signals
- `simulation_plots.m`
    - Plots bar graphs and creates tables for percentage of algorithm runs with ratio of 2nd hughest to highest peak frequency amplitudes after separation


## Array Coherence Folder

- `matrix_reordering.m`
    - Plots sparsity structure of matrix with SRCM and SAMD permutations
- `coherence_square_hat.m`
    - Calculates coherence matrix for square and hat function
    - Performs SRCM and SAMD matrix permutations
- `coherence_square_hat_shift.m`
    - Calculates coherence matrix for shifted square and hat functions
    - Performs SRCM and SAMD matrix permutations
- `coherence_square_hat_shift_nonideal.m`
    - Calculates coherence matrix for shifted square and hat function with nonideal separation
    - Performs SRCM and SAMD matrix permutations
- `coherence_real_data.m`
    - Calculates coherence matrix for real data loaded in 
    - Performs SRCM and SAMD matrix permutations
    
    
## Test Data Folder

- where to put data to be loaded in for testing of real data examples



## k-Wave

This code uses V1.3 of the Matlab k-Wave toolbox for a couple of tests in the ICA folder.


B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave-fields," J. Biomed. Opt., vol. 15, no. 2, p. 021314, 2010.

