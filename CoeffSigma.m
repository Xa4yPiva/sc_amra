function [sigmaAP, sigmaDP] = CoeffSigma(signalHilb, aNorm, aThreshold, fc, fs)
%COEFFSIGMA Summary of this function goes here
%   Detailed explanation goes here

phi = angle(signalHilb);
phiUW = unwrap(phi);
phiNoLin = detrend(phiUW);
i = (0 : length(phiNoLin)-1);
phiNoLin1 = phiUW - (2*pi*fc*i)/fs;

phiNL = phiNoLin(aNorm > aThreshold);
phiNL1 = phiNoLin1(aNorm > aThreshold);
C = length(phiNL);
sumPhiNL2 = sum(phiNL .* phiNL);
sumAbsPhiNL = sum(abs(phiNL));
sumPhiNL = sum(phiNL);
sigmaAP = sqrt((sumPhiNL2/C) - ((sumAbsPhiNL/C) * (sumAbsPhiNL/C)))
sigmaDP = sqrt((sumPhiNL2/C) - ((sumPhiNL/C) * (sumPhiNL/C)))

sumPhiNL12 = sum(phiNL1 .* phiNL1);
sumAbsPhiNL1 = sum(abs(phiNL1));
sumPhiNL1 = sum(phiNL1);
sigmaAP1 = sqrt((sumPhiNL12/C) - ((sumAbsPhiNL1/C) * (sumAbsPhiNL1/C)))
sigmaDP1 = sqrt((sumPhiNL12/C) - ((sumPhiNL1/C) * (sumPhiNL1/C)))

% figure(3);
% subplot(4,1,1); plot(phiUW);        grid on;
% subplot(4,1,2); plot(phiNoLin);     grid on;
% subplot(4,1,3); plot(phiNL);        grid on;
% subplot(4,1,4); plot(abs(phiNL));   grid on;

% figure(4);
% subplot(3,1,1); plot(phiUW);        grid on;
% subplot(3,1,2); plot(phiNoLin);     grid on;
% subplot(3,1,3); plot(phiNoLin1);     grid on;

end

