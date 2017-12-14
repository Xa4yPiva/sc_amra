function [sigmaAP, sigmaDP] = CoeffSigma(envelope, aNorm, aThreshold)
%COEFFSIGMA Summary of this function goes here
%   Detailed explanation goes here


% phi = angle(envelope*exp(1i*(-pi/2)));
% % phiUW = unwrap(phi);
% phiUW = phi;
% phiCent = phiUW - mean(phiUW);
% % phiCentAbs = abs(phiUW) - mean(abs(phiUW));
% phiNL = phiCent(aNorm > aThreshold);
% phiN = phiNL;

phi = angle(envelope*exp(1i*(-pi/2)));
phiNW = phi(aNorm > aThreshold);
% phiUW = unwrap(phiNW);
phiUW = phiNW;
phiCent = phiUW - mean(phiUW);
phiN = phiCent;

% figure(3);
% subplot(2,2,1);
% plot(angle(envelope));
% grid on;
% subplot(2,2,2);
% plot(phi);
% grid on;
% subplot(2,2,3);
% plot(phiUW);
% grid on;
% subplot(2,2,4);
% plot(phiCent);
% grid on;

C = length(phiN);
sumPhiNL2 = sum(phiN .* phiN);
sumAbsPhiNL = sum(abs(phiN));
sumPhiNL = sum(phiN);
sigmaAP = sqrt((sumPhiNL2/C) - ((sumAbsPhiNL/C) * (sumAbsPhiNL/C)));
sigmaDP = sqrt((sumPhiNL2/C) - ((sumPhiNL/C) * (sumPhiNL/C)));

end

