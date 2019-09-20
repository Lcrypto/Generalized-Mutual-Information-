function [GMI,MIperBitPosition]=calcGMI(X,Y,labeling)
%calcMI calculates the sum of bitwise memoryless mutual informations using circularly symmetric Gaussian noise statistics.
% This quantity is also known as generalized mutual information (GMI). In
% contrast to the symbolwise MI, it is an achievable rate for receivers with binary decoding and no iterations between demapper and decoder.
% The current version of this script works only for uniformly distributed square QAM.
%
% Input:
% X         1 x N       N transmitted complex symbols chosen from 2^m possible values
% Y         1 x N       N received complex symbols
% labeling  chars       optional: 'Gray' (default) or 'Binary'
%
% Output:
% GMI                   1 x 1       Memoryless generalized mutual information assuming circularly symmetric Gaussian noise statistics
% MIperBitPosition      m x 1       Mutual information of each of the m parallel bit channels
%
% Author: Tobias Fehenberger <tobias.fehenberger@tum.de>, Apr. 2015

%% number of bits per symbol
m = log2(numel(unique(X)));

%% we want Y as column vector and X as row vector
if size(Y,1) ~= 1
    Y = Y.';
end
if size(X,2) ~= 1
    X = X.';
end

%% we need X also in integer representation
% The input is assumed to be square QAM. For other formats, adapt the demodulator object.
X = X/sqrt(mean(abs(X).^2)); % normalize
hDemod = comm.RectangularQAMDemodulator(2^m, 'BitOutput',false, ...
    'NormalizationMethod', 'Average power', 'SymbolMapping', 'Binary'); % mapping is only important for shaped X when
                                                                        % the symbols must be sorted from top left
                                                                        % running down column-wise to the bottom right.
Xint = step(hDemod,X)';

%% Labeling
if nargin==2
    labeling = 'Gray';
end

%% Calculate sent bits and LLRs
Y = Y/sqrt(mean(abs(Y).^2));
N0 = estimateNoiseVar(Xint,Y);

hDemodInput = comm.RectangularQAMDemodulator(2^m, 'BitOutput',true, ...
    'NormalizationMethod', 'Average power', ...
    'SymbolMapping', labeling);

hDemodOutput = comm.RectangularQAMDemodulator(2^m, 'BitOutput',true, ...
    'DecisionMethod', 'Log-likelihood ratio', 'SymbolMapping', labeling,...
    'Variance', N0, 'NormalizationMethod', 'Average power');

c = step(hDemodInput, X)';
LLRs = step(hDemodOutput,Y.')';

%% Compute bitwise MIs and their sum
MIperBitPosition=zeros(m,1);

for kk=1:m
    MIperBitPosition(kk)=1-mean(log2(1+exp((2.*c(kk:m:end)-1).*LLRs(kk:m:end))));
end
GMI=sum(MIperBitPosition);
end

function [N0,mus,inputPower,SNR,P_X] = estimateNoiseVar(Xint,Y)
%estimateNoiseVar estimates the 2D (complex) noise variance of N received symbols
% The estimate comes from Matlab's built-in normfit function giving the mean and the standard deviation (NOT the variance!).
%
% Input:
% Xint          1 x N         sent symbols in integer representation from 0 to M-1
% Y             1 x N         complex received symbols. Normalized.
%
% Output:
% N0            1 x 1       Twice the noise variance in each dimension
% mus           M x 1       Complex mean of each constellation cloud
% inputpower    1x1       	Input power
% SNR           1 x 1    	SNR in dB: 10*log10(inputPower/N0)
% P_X           M x 1       Probability mass function of input
%
% Author: Tobias Fehenberger <tobias.fehenberger@tum.de>, Jan. 2015

M = numel(unique(Xint));

N0 = 0;
inputPower = 0;
mus = zeros(M,1);
P_X = zeros(M,1);

for m=0:M-1
    mthSymbol = Y(Xint==m);
    P_X(m+1) = length(mthSymbol)/length(Y);
    [mu,mthsigmahat] = normfit(mthSymbol);
    mus(m+1) = mu;
    N0 = N0+(mthsigmahat)^2*P_X(m+1);
    inputPower = inputPower + abs(mu)^2*P_X(m+1);
end
SNR = 10*log10(inputPower/N0);
end