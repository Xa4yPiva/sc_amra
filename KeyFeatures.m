function [P, gammaMax, sigmaAP, sigmaDP] = KeyFeatures(signal, fc, fs, aThreshold)
%KEYFEATURES Summary of this function goes here
%   Detailed explanation goes here

len = length(signal);
Nfft = 2 ^ nextpow2(len);

z = hilbert(signal, Nfft);

P = CoeffP(signal, Nfft, fc, fs);
[aNorm, aCN] = SubCoeffsA(z);
gammaMax = CoeffGammaMax(aCN, Nfft);
[sigmaAP, sigmaDP] = CoeffSigma(z, aNorm, aThreshold, fc, fs);


end

