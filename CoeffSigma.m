function [sigmaAP, sigmaDP] = CoeffSigma(envelope, aNorm, aThreshold)
%COEFFSIGMA Summary of this function goes here
%   Detailed explanation goes here

phi = angle(envelope);
phiNL = phi(aNorm > aThreshold);
C = length(phiNL);
sumPhiNL2 = sum(phiNL .* phiNL);
sumAbsPhiNL = sum(abs(phiNL));
sumPhiNL = sum(phiNL);
sigmaAP = sqrt((sumPhiNL2/C) - ((sumAbsPhiNL/C) * (sumAbsPhiNL/C)));
sigmaDP = sqrt((sumPhiNL2/C) - ((sumPhiNL/C) * (sumPhiNL/C)));

end

