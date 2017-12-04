clear all;

fs = 1200e3;
fc = 150e3;
f1 = 1000;
f2 = 3000;
n = 2048;
dt = 1 / fs;
t = (0 : n-1) * dt;
N = 2^(nextpow2(n)+2);
f = (-N/2 : (N-1)/2) * fs / N;
x = 1*cos(2 * pi * f1 * t);% + 0.5*cos(2 * pi * f2 * t);

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
% subplot(1,2,1); plot(t, xFM);       grid on;
% subplot(1,2,2); plot(f, specFM);    grid on;

offset = exp(1i * 2*pi*(-fc) * t);

signal = awgn(xAM, 20);
% signal = xAM;

specSig = fftshift(abs(fft(signal, N))) / (N/2);
envelope = hilbert(signal) .* offset;
specEnv = fftshift(abs(fft(envelope, N))) / (N/2);

figure(1)
subplot(2,2,1); plot(t, signal);       grid on;
subplot(2,2,3); plot(f, specSig);    grid on;
subplot(2,2,2); plot(t, envelope);       grid on;
subplot(2,2,4); plot(f, specEnv);   grid on;

phi = angle(hilbert(signal));
phiUW = unwrap(phi);
phiNoLin = detrend(phiUW);
phi1 = angle(envelope);
phiUW1 = unwrap(phi1);
phiNoLin1 = detrend(phiUW1);

figure(5);
subplot(2,2,1); plot(phiUW);        grid on;
subplot(2,2,2); plot(phiNoLin);     grid on;
subplot(2,2,3); plot(phiUW1);        grid on;
subplot(2,2,4); plot(phiNoLin1);     grid on;

figure(6);
plot(angle(envelope)); grid on;

%----------------AMRA-------------------
atopt = 1; % super dirty hardcode! Compute threshold or something!
[P, yMax, sAP, sDP] = KeyFeatures(signal, fc, fs, atopt);

% Super dirty hardcode!
threshP = 0.5;
threshY = 4;
threshSap = 2*pi;
threshSdp = 0.15;

amra1 = AMRA1(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);
amra2 = AMRA2(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);
amra3 = AMRA3(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);
amra4 = AMRA4(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);
amra5 = AMRA5(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp);

fprintf("key features:\nP = %f\nyMax = %f\nsAP = %f\nsDP = %f\n", P, yMax, sAP, sDP);
fprintf("amra1 = " + amra1 + "\n");
fprintf("amra2 = " + amra2 + "\n");
fprintf("amra3 = " + amra3 + "\n");
fprintf("amra4 = " + amra4 + "\n");
fprintf("amra5 = " + amra5 + "\n");


% figure(2);
% subplot(2,1,1); plot(signal);                 grid on;
% subplot(2,1,2); plot(abs(hilbert(signal)));   grid on;






