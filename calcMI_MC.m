function [MI,N0,SNR]=calcMI_MC(X,Y,N0)
%calcMI_MC calculates the symbolwise memoryless mutual information (MI) using circularly symmetric Gaussian noise statistics.
% The estimate is calculated in a Monte Carlo approach using n symbols.
% If the channel is not circularly symmetric Gaussian, a lower bound on MI, i.e., an achievable rate, is obtained.
%
% To have a reliable estimate, you should first obtain the noise variance
% N0 from a first set of samples. In a second step, estimate the MI on a
% different sequence of symbols with the previously found noise variance.
% This corresponds to a double Monte Carlo approach.
%
% The $M$ means of the received constellation points are assumed to
% correspond to the input $X$. Adaptive partitioning is thus not employed.

% Input:
% X     1 x N       N transmitted complex symbols, uniformly distributed
% Y     1 x N       N received complex symbols
% N0    1 x 1       optional: 2D noise variance for double Monte Carlo

%
% Output:
% MI    1 x 1       Memoryless mutual information assuming circularly symmetric Gaussian noise statistics
% N0    1x1         noise variance
%
% Author: Tobias Fehenberger <tobias.fehenberger@tum.de>, Aug. 2015

%% Input Manipulation
% we want X and Y as row vectors
if size(Y,1) ~= 1
    Y = Y.';
end
if size(X,1) ~= 1
    X = X.';
end

%% Variable Initialization
% find modulation order M and initialize variables
M = numel(unique(X));
if ~~(fix(log2(M))-log2(M)) %#ok<BDLOG> % if M is not a power of two
    M=2^ceil(log2(M));
end
n=length(X);
constellation=zeros(M,1);
P_X=zeros(M,1);

%% Get X in Integer Representation
X = X/sqrt(mean(abs(X).^2)); % normalize such that var(X)=1
Y = Y/sqrt(mean(abs(Y).^2)); % normalize such that var(Y)=1
avgPower=(min(unique(abs(X)))/abs(1/sqrt(2/3*(M-1))+1j*1/sqrt(2/3*(M-1))))^2;
% uniform input: avgPower equals 1
% shaped input: avgPower is the increase in minium squared Euclidean distance
hDemod = comm.RectangularQAMDemodulator(M, 'BitOutput',false, ...
    'NormalizationMethod', 'Average power', 'AveragePower', avgPower, ...
    'SymbolMapping', 'Binary'); % mapping is only important for shaped X when
% the symbols must be sorted from top left
% running down column-wise to the bottom right.
Xint = step(hDemod,X.')';

%% If the noise variance is not given, calculate it and scale X such that the means of X and Y are aligned
if nargin==2
    fun = @(h) (h*X-Y)*(h*X-Y)';
    h=fminbnd(fun,0,2);
    N0=(1-h^2)/h^2;
else
    h=sqrt(1/(N0+1));
end
Y=Y/h;

%% Find constellation and empirical input distribution
for s=0:M-1
    constellation(s+1,:)=X(:,find(Xint==s,1));
    P_X(s+1)=nnz(Xint==s)/n;
end

%% Monte Carlo estimation of (a lower bound to) the mutual information I(X;Y)
qYonX=(1/(pi*N0)*exp((-(real(Y)-real(X)).^2-(imag(Y)-imag(X)).^2)/N0));
qY=0;
for ii=1:M
    qY=qY+P_X(ii)*(1/(pi*N0)*exp((-(real(Y)-real(constellation(ii))).^2-(imag(Y)-imag(constellation(ii))).^2)/N0));
end
MI=1/n*sum(log2(max(qYonX,realmin)./max(qY,realmin)));

end