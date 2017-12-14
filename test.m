clear all;

fs = 1200e3;
fc = 150e3;
f1 = 1000;
f11 = f1 * 2;
f2 = 8000;
n = 2048;
dt = 1 / fs;
t = (0 : n-1) * dt;
N = 2^(nextpow2(n)+4);
f = (-N/2 : (N-1)/2) * fs / N;
% x = 0.5*cos(2 * pi * f1 * t);
x = 0.5*cos(2 * pi * f1 * t) + 0.25*cos(2 * pi * (f11) * t);

xAM = ammod(x, fc, fs, 0, 1);
specAM = fftshift(abs(fft(xAM, N))) / (N/2);
xDSB = ammod(x, fc, fs, 0, 0);
specDSB = fftshift(abs(fft(xDSB, N))) / (N/2);
xLSB = ssbmod(x, fc, fs, 0);
specLSB = fftshift(abs(fft(xLSB, N))) / (N/2);
xUSB = ssbmod(x, fc, fs, 0, 'upper');
specUSB = fftshift(abs(fft(xUSB, N))) / (N/2);
fDev = f2 - f1;
xFM = fmmod(x, fc, fs, fDev);
specFM = fftshift(abs(fft(xFM, N))) / (N/2);
mComb = 0.5;
xComb = (1 + mComb*x) .* xFM;
specComb = fftshift(abs(fft(xComb, N))) / (N/2);

% figure(1);
% subplot(2,2,1); plot(t, xAM);       grid on;
% subplot(2,2,2); plot(f, specAM);    grid on;
% subplot(2,2,3); plot(t, xDSB);      grid on;
% subplot(2,2,4); plot(f, specDSB);   grid on;
% figure(2);
% subplot(2,2,1); plot(t, xLSB);      grid on;
% subplot(2,2,2); plot(f, specLSB);   grid on;
% subplot(2,2,3); plot(t, xUSB);      grid on;
% subplot(2,2,4); plot(f, specUSB);   grid on;
% figure(3);
% subplot(2,2,1); plot(t, xFM);       grid on;
% subplot(2,2,2); plot(f, specFM);    grid on;
% subplot(2,2,1); plot(t, xComb);     grid on;
% subplot(2,2,2); plot(f, specComb);  grid on;

offset = exp(1i * 2*pi*(-fc) * t);

signal = awgn(xDSB, 10);
% signal = xDSB;
envel = hilbert(signal) .* offset;
specSig = fftshift(abs(fft(signal, N))) / (N/2);
specEnv = fftshift(abs(fft(envel, N))) / (N/2);

% figure(2);
% subplot(3,1,1);
% plot(signal);
% grid on;
% subplot(3,1,2);
% plot(abs(envel));
% grid on;
% subplot(3,1,3);
% plot(angle(envel*1*exp(1i*(-pi/2))));
% % plot(angle(envel(2:end) .* conj(envel(1:end-1))));
% grid on;

%----------------AMRA-------------------
atopt = 0.5; % dirty hardcode! Compute threshold or something!
[P, yMax, sAP, sDP] = KeyFeatures(envel, atopt);

% Dirty hardcode!
threshP = 0.5;
threshY = 5.5;
threshSap = pi/4;
threshSdp = pi/6;

amra1 = AMRA1(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);
% amra2 = AMRA2(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);
% amra3 = AMRA3(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);
% amra4 = AMRA4(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);
% amra5 = AMRA5(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);

fprintf("key features:\nP = %f\nyMax = %f\nsAP = %f\nsDP = %f\n", P, yMax, sAP, sDP);
fprintf("amra1 = " + amra1 + "\n");


%---------Probability-of-right-decision-----------
signals = [xAM; xDSB; xLSB; xUSB; xFM; xComb];
decRight = ["AM", "DSB", "LSB", "USB", "FM", "Combined"];
sigsNum = length(decRight);
SNR = (-6 : 0.5 : 8);
lenSNR = length(SNR);
expNum = 1000;
pRight = zeros(sigsNum, lenSNR);
cyclesNum = sigsNum*lenSNR;
iteration = 0;
h = waitbar(0, 'Computing...');
tic
for k = 1 : sigsNum
    envk = hilbert(signals(k,:)) .* offset;
    for i = 1 : lenSNR
        amra1_1 = zeros(1, expNum);
        for j = 1 : expNum
            env = awgn(envk, SNR(i));
            [P1, yMax1, sAP1, sDP1] = KeyFeatures(env, atopt);
            decision = AMRA1(P1, yMax1, sAP1, sDP1, threshP, threshY, threshSap, threshSdp);
            if (decision == decRight(k))
                amra1_1(j) = 1;
            end
        end
        pRight(k, i) = sum(amra1_1) / expNum;
        iteration = iteration + 1;
        waitbar(iteration / cyclesNum);
    end
end
toc
close(h);

mType = ['x', '*', 's', '^', 'o', 'd-'];
figure(1);
for i = 1 : sigsNum
    plot(SNR, pRight(i,:), 'marker', mType(i), 'markersize', 10, 'linewidth', 2);
    hold on;
end
title('AMRA1')
xlabel('SNR, dB');
ylabel('Probability of right decision');
grid on;
legend(decRight, 'location', 'northwest');
legend('show');
set(groot, 'DefaultAxesFontSize', 18);


