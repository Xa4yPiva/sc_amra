function [aNorm, aCentNorm] = SubCoeffsA(envelope)
%SUBCOEFFSA Summary of this function goes here
%   Detailed explanation goes here

a = abs(envelope);
aMean = mean(a);
aNorm = a / aMean;
aCentNorm = aNorm - 1;
% figure(2); plot(aCentNorm); grid on;

end

