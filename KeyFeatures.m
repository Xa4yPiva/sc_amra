function [P, gammaMax, sigmaAP, sigmaDP] = KeyFeatures(envelope, aThreshold)
%KEYFEATURES Summary of this function goes here
%   Detailed explanation goes here

len = length(envelope);
Nfft = 2 ^ nextpow2(len);

P = CoeffP(envelope, Nfft);
[aNorm, aCN] = SubCoeffsA(envelope);
gammaMax = CoeffGammaMax(aCN, Nfft);
[sigmaAP, sigmaDP] = CoeffSigma(envelope, aNorm, aThreshold);

end

