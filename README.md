# Generalized-Mutual-Information
The GitHub repository contains tools for estimating achievable rates, including Mutual Information (MI) and Generalized Mutual Information (GMI), developed by Dr. Tobias Fehenberger (https://www.fehenberger.de).

When deciding which metric to use, if you have a given input and a memoryless channel, MI is the largest achievable rate. However, if you want to restrict the receiver to binary decoding and iterations between demapper and decoder are not allowed, then GMI becomes an achievable rate (MI is generally not).

The repository also includes a tool for calculating Signal-to-Noise Ratio (SNR) and Optical Signal-to-Noise Ratio (OSNR) based on the GN model for multispan Nyquist-Wavelength Division Multiplexing (WDM) with Erbium-Doped Fiber Amplifiers (EDFAs). This tool implements Equations (16), (22), (36), (50), and (40) from P. Poggiolini, G. Bosco, and A. Carena's paper titled "The GN-Model of Fiber Non-Linear Propagation and its Applications," published in the Journal of Lightwave Technology in February 2014.

