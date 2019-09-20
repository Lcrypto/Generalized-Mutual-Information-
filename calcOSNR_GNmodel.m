% Calculates SNR and OSNR according to the GN model for multispan Nyquist-WDM with EDFAs. Implementated are Eqs. (16), (22), (36), (50), and (40) of [1]
% The results are stored in the struct GN. 
%
% Author: Tobias Fehenberger <tobias.fehenberger@tum.de>, Apr. 2015
%
% [1] P. Poggiolini, G. Bosco, and A. Carena, "The GN-Model of Fiber Non-Linear Propagation and its Applications," J. Light. Technol., vol. 32, no. 4, pp. 694-721, Feb. 2014.

%% System parameters
D=17;         % ps/nm/km
gamma=1.3e-3;	% 1/W/m
B_signal=30e9;	% signal bandwidth = channel spacing for Nyquist WDM [Hz]
N_channel=17;	% numbers of channels
NFdB=4;         % amplifier noise figure [dB]
Length=100e3;	% span length (m)
alpha=.2e-3;	% fibre attenuation (dB/m)
N_spans=10;     % number of spans
N_pol=2;        % number of polarizations

%% Constants
h=6.6256e-34;	% Planck's constant [J/s]
nu=299792458/1550e-9;	% light frequency [Hz]
SNR2OSNR = 10*log10(B_signal/12.5e9*N_pol/2); % Eq. (34) of Essiambre et al., "Capacity Limits of Optical Fiber Networks," J. Light. Technol., vol. 28, no. 4, pp. 662-701, Feb. 2010.
dB2Neper = 10/log(10);

%% Set launch power per channel
powerVec = -10:.25:10; % launch power [dBm]

%% Some more quantities
B_total = N_channel*B_signal; % total system bandwidth
beta2 = -(1550e-9)^2 * (D*1e-6) / (2*pi*3e8);	% propagation constant
GaindB = Length*alpha;	% amplifier gain (dB)
L_eff = ((1-exp(-(alpha/dB2Neper)*Length))/(alpha/dB2Neper)); % effective length [m]
L_effa = 1/(alpha/dB2Neper); % asymptotic effective length [m]
G_tx_ch = 10.^((powerVec-30)/10)./(B_signal); % [W/Hz]
ASE = N_pol*N_spans*10^(NFdB/10)/2*(10^(GaindB/10)-1)*h*nu; % Eq. (50)


%% GN model
epsilon = 3/10*log(1+ 6/Length * L_effa / (asinh(0.5*pi^2*abs(beta2)*L_effa*B_total^2))); % Eq. (40)
G_NLI = gamma^2.*G_tx_ch.^3*L_eff^2*(2/3)^3*asinh(0.5*pi^2*abs(beta2)*L_effa*B_total^2)/(pi*abs(beta2)*L_effa); % Eq. (36)

%% SNR calculation
GN.SNR_NLI = 10*log10(G_tx_ch ./ (ASE + N_spans^(1+epsilon)*G_NLI)); % Eq. (16), (22)
GN.SNR_ASE = 10*log10(G_tx_ch ./ ASE);

%% OSNR calculation
GN.OSNR_NLI = GN.SNR_NLI+SNR2OSNR;
GN.OSNR_EDFA = GN.SNR_ASE+SNR2OSNR;

GN.power = powerVec;