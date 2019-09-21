# Generalized-Mutual-Information
Achievable Rates: Mutial Information (MI) and Generalize Mutial Information (GMI)
from Dr. Tobias Fehenberger https://www.fehenberger.de


When to use (symbolwise) MI and when to use GMI


For a given input and a memoryless channel, the MI is the largest achievable rate. If you do not want to make restrictions on the receiver and its decoder, use MI. If you consider a receiver with binary decoding and iterations between demapper and decoder are not allowed, GMI is an achievable rate (MI is in general not!).


Calculates SNR and OSNR according to the GN model for multispan Nyquist-WDM with EDFAs. Implementated are Eqs. (16), (22), (36), (50), and (40) of  P. Poggiolini, G. Bosco, and A. Carena, "The GN-Model of Fiber Non-Linear Propagation and its Applications," J. Light. Technol., vol. 32, no. 4, pp. 694-721, Feb. 2014.

