addpath(genpath('../common'));
clear;
close all;


%% Signal
fs0 = 12e3;
T = 0.1;
fLowHz = 300;
fHighHz = 3600;
bandHz = fHighHz - fLowHz;
data = RandomBandLimitedSignal(fs0, T, 20, fLowHz, fHighHz, 4000, 60, 1, 60, 'uniform');
data1 = RandomBandLimitedSignal(fs0, T, 20, fLowHz, fHighHz/2, 4000/2, 60, 1, 60, 'uniform');
fs = 1200e3;
factor = fs / fs0;
[p, q] = rat(factor);
data = resample(data, p, q);
data1 = resample(data1, p, q);
lenFrame = 8192;
framesNum = floor(length(data) / lenFrame);
x = data;
x1 = data1;
lenSignal = length(x);
N = 2 ^ nextpow2(lenSignal);
fc = 60e3;
dt = 1 / fs;
t = (0 : lenSignal-1) * dt;
f = (-N/2 : (N-1)/2) * fs / N;
specX = fftshift(abs(fft(x, N))) / (N/2);

offset = exp(1i * 2*pi*fc * t);

%% Modulation
mAM = 0.5;
xAM = ammod(mAM * x, fc, fs, 0, 0.5);
specAM = fftshift(abs(fft(xAM, N))) / (N/2);
xDSB = ammod(x, fc, fs, 0, 0);
specDSB = fftshift(abs(fft(xDSB, N))) / (N/2);
xLSB = ssbmod(x, fc, fs, 0);
specLSB = fftshift(abs(fft(xLSB, N))) / (N/2);
xUSB = ssbmod(x, fc, fs, 0, 'upper');
specUSB = fftshift(abs(fft(xUSB, N))) / (N/2);
% fDev = f2 - f1;
fDev = 25e3;
xFM = fmmod(x, fc, fs, fDev);
specFM = fftshift(abs(fft(xFM, N))) / (N/2);
mComb = 0.5;
xComb = (1 + mComb*x1) .* xFM;
specComb = fftshift(abs(fft(xComb, N))) / (N/2);

uvsbFilterSpec = fdesign.highpass('N,Fc,Ast,Ap', 1000, fc, 60, 1, fs);
uvsbFilter = design(uvsbFilterSpec, 'equiripple');
lvsbFilterSpec = fdesign.lowpass('N,Fc,Ap,Ast', 1000, fc, 1, 60, fs);
lvsbFilter = design(lvsbFilterSpec, 'equiripple');
hu = freqz(uvsbFilter, length(f)/2);
hl = freqz(lvsbFilter, length(f)/2);
% figure(2); set(gcf, 'color', 'w');
% plot(f(end/2+1 : end), abs(hu))
% % hold on; plot(f(end/2+1 : end), abs(hl));
% grid on;
% xlabel('f, Hz'); ylabel('|H(f)|');
xUVSB = filter(uvsbFilter, xAM);
xLVSB = filter(lvsbFilter, xAM);
specUVSB = fftshift(abs(fft(xUVSB, N))) / (N/2);
specLVSB = fftshift(abs(fft(xLVSB, N))) / (N/2);

%% Plot signals
% figure(1); set(gcf, 'color', 'w');
% subplot(2,1,1); plot(t, x);             grid on; xlabel('t, s'); ylabel('s(t)'); title('Noise in speech band');
% subplot(2,1,2); plot(f, mag2db(specX)); grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([-5000, 5000]);

% figure(2); set(gcf, 'color', 'w');
% subplot(2,2,1); plot(t, xAM);               grid on;
% subplot(2,2,2); plot(f, mag2db(specAM));    grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]);
% subplot(2,2,3); plot(t, xDSB);              grid on;
% subplot(2,2,4); plot(f, mag2db(specDSB));   grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]);
% figure(3); set(gcf, 'color', 'w');
% subplot(2,2,1); plot(t, xLSB);              grid on;
% subplot(2,2,2); plot(f, mag2db(specLSB));   grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]);
% subplot(2,2,3); plot(t, xUSB);              grid on;
% subplot(2,2,4); plot(f, mag2db(specUSB));   grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]);
% figure(4); set(gcf, 'color', 'w');
% subplot(2,2,1); plot(t, xFM);               grid on;
% subplot(2,2,2); plot(f, mag2db(specFM));    grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-fDev, fc+fDev]);
% subplot(2,2,3); plot(t, xComb);             grid on;
% subplot(2,2,4); plot(f, mag2db(specComb));  grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-fDev, fc+fDev]);
% figure(5); set(gcf, 'color', 'w');
% subplot(2,2,1); plot(t, xUVSB);             grid on;
% subplot(2,2,2); plot(f, mag2db(specUVSB));  grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]);
% subplot(2,2,3); plot(t, xLVSB);             grid on;
% subplot(2,2,4); plot(f, mag2db(specLVSB));  grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]);

% figure(2); set(gcf, 'color', 'w');
% subplot(2,1,1); plot(f, mag2db(specAM));  grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]); title('AM');
% subplot(2,1,2); plot(f, mag2db(specDSB)); grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]); title('DSB');
% figure(3); set(gcf, 'color', 'w');
% subplot(2,1,1); plot(f, mag2db(specLSB)); grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]); title('LSB');
% subplot(2,1,2); plot(f, mag2db(specUSB)); grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]); title('USB');
% figure(4); set(gcf, 'color', 'w');
% subplot(2,1,1); plot(f, mag2db(specFM));   grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-2*fDev, fc+2*fDev]); title('FM');
% subplot(2,1,2); plot(f, mag2db(specComb)); grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-2*fDev, fc+2*fDev]); title('Combined (FM + AM)');
% figure(5); set(gcf, 'color', 'w');
% subplot(2,1,1); plot(f, mag2db(specUVSB)); grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]); title('LVSB');
% subplot(2,1,2); plot(f, mag2db(specLVSB)); grid on; xlabel('f, Hz'); ylabel('|S(f)|, dB'); xlim([fc-5000, fc+5000]); title('UVSB');

%% Envelope
% signal = awgn(xDSB, 10);
signal = xFM(328 : 328 + lenFrame-1);
envel = hilbert(signal) .* conj(offset(1 : lenFrame));
envel = awgn(envel, 10, 'measured');
specSig = fftshift(abs(fft(signal, N))) / (N/2);
specEnv = fftshift(abs(fft(envel, N))) / (N/2);
% figure(5);
% subplot(2,1,1); plot(f, mag2db(specSig)); grid on;
% subplot(2,1,2); plot(f, mag2db(specEnv)); grid on;

%% AMRA OneShotTest
% Dirty hardcode! Compute threshold or something!
thresholds.ampl = 1;
thresholds.P = 0.5;
thresholds.gammaMax = 5.5;
thresholds.sigmaAP = pi/4;
thresholds.sigmaDP = pi/3;

kf = KeyFeatures(envel, thresholds.ampl);

amra1 = AMRA1(kf, thresholds);
% amra2 = AMRA2(kf, thresholds);
% amra3 = AMRA3(kf, thresholds);
% amra4 = AMRA4(kf, thresholds);
% amra5 = AMRA5(kf, thresholds);

fprintf("key features:\nP = %f\nyMax = %f\nsAP = %f\nsDP = %f\n", kf.P, kf.gammaMax, kf.sigmaAP, kf.sigmaDP);
fprintf("amra1 = " + amra1 + "\n");
% fprintf("amra2 = " + amra2 + "\n");
% fprintf("amra3 = " + amra3 + "\n");
% fprintf("amra4 = " + amra4 + "\n");
% fprintf("amra5 = " + amra5 + "\n");

%% Signals, Decisions
bandSignalMod = 3300;
signals = [xAM; xDSB; xLSB; xUSB; xLVSB; xUVSB; xFM; xComb];
bandsHz = [2*fHighHz, 2*fHighHz, bandHz, bandHz, 1.2*bandHz, 1.2*bandHz, 2*fDev, 2*fDev];
sigsNum = min(size(signals));
decRight = ["AM", "DSB", "LSB", "USB", "LVSB", "UVSB", "FM", "Combined"];
mType = ['x', '*', 's', '^', 'p', 'h', 'o', 'd'];

%% Probability of right decision
envelopes = zeros(sigsNum, lenSignal);
energies = zeros(1, sigsNum);
for i = 1 : sigsNum
    envelopes(i, :) = hilbert(signals(i, :)) .* conj(offset);
%     energies(i) = sum(abs(envelopes(i, :)) .^ 2) / fs;
%     envelopes(i, :) = envelopes(i, :) / sqrt(energies(i));
end
SNR = (-8 : 1 : 10);
expNum = 100;
pRight = ProbRightDecision(envelopes, thresholds, decRight, SNR, expNum, lenFrame);

figure(6);
for i = 1 : sigsNum
    plot(SNR, pRight(i,:), 'marker', mType(i), 'markersize', 10, 'linewidth', 2);
    hold on;
end
grid on;
title('AMRA1'); xlabel('SNR, dB'); ylabel('Probability of right decision');
legend(decRight, 'location', 'northwest'); legend('show');
set(gcf, 'color', 'w'); set(groot, 'DefaultAxesFontSize', 18);

%% Offset Test
% modulation = "Combined";
% idxRight = find(decRight == modulation);
% 
% % fOff = -bandHz : bandHz;
% % fOff = -5*fs/N : fs/N : 5*fs/N;
% % fOff = linspace(-bandHz*2, bandHz*2, 21);
% % fOff = linspace(-bandHz/2, bandHz/2, 11);
% % fOff = linspace(-100, 100, 21);
% fOff = linspace(-bandsHz(idxRight)/4, bandsHz(idxRight)/4, 41);
% lenOff = length(fOff);
% fOffset = zeros(lenOff, lenSignal);
% for i = 1 : lenOff
%     fOffset(i, :) = conj(offset) .* exp(1i * 2*pi*(-fOff(i)) * t);
% end
% 
% env1 = hilbert(signals(idxRight, :));
% 
% snr = 10;
% expNum = 100;
% cyclesNum = lenOff * expNum;
% 
% prob=zeros(sigsNum, lenOff);
% iteration = 0;
% h = waitbar(0, 'Computing...');
% for i = 1 : lenOff
%     envi = env1 .* fOffset(i,:);
%     for j = 1 : expNum
% %         pos = mod(j, framesNum);
%         pos = floor(rand() * (framesNum-1));
%         env = envi(pos*lenFrame+1 : (pos+1)*lenFrame);
%         env = awgn(env, snr, 'measured');
%         kf = KeyFeatures(env, thresholds.ampl);
%         decision = AMRA1(kf, thresholds);
%         idx = find(decRight == decision);
%         prob(idx, i) = prob(idx, i) + 1;
%         iteration = iteration + 1;
%         waitbar(iteration / cyclesNum);
%     end
% end
% close(h);
% 
% prob = prob / expNum;
% 
% figure(7);
% for i = 1 : sigsNum
%     if i == idxRight
%         line = '-';
%     else
%         line = '--';
%     end
% %     plot(fOff, prob(i,:), line, 'marker', mType(i), 'markersize', 10, 'linewidth', 2);
%     plot(fOff/bandsHz(idxRight), prob(i,:), line, 'marker', mType(i), 'markersize', 10, 'linewidth', 2);
%     hold on;
% end
% grid on;
% title(strcat(modulation, ", ", num2str(snr), " dB"));
% ylabel('Probability of right decision');
% % xlabel('\Delta f, Hz');
% xlabel('$\displaystyle\frac{\delta f}{F}$','interpreter','latex');
% legend(decRight); legend('show');
% set(groot, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');

%% Weights
modulation = "AM";
idxRight = find(decRight == modulation);
envelope = hilbert(signals(idxRight, :)) .* conj(offset);

snr = [-2.5, -5, 0.5, 0.5, 5, 5, 0, 2];
lenFrame = 8192;
% framesNum = floor(lenSignal / lenFrame);
% w = zeros(1, sigsNum);
% cyclesNum = framesNum;
% iteration = 0;
% h = waitbar(0, 'Computing...');
% for i = 1 : framesNum
%     env = envelope((i-1)*lenFrame+1 : i*lenFrame);
%     env = awgn(env, snr, 'measured');
%     kf = KeyFeatures(env, thresholds.ampl);
%     decision = AMRA1(kf, thresholds);
%     idx = find(decRight == decision);
%     w(idx) = w(idx) + 1;
%     iteration = iteration + 1;
%     waitbar(iteration / cyclesNum);
% end
% close(h);
% weights = w / framesNum;

weights = zeros(sigsNum, sigsNum);
for i = 1 : sigsNum
    envelope = hilbert(signals(i, :)) .* conj(offset);
    weights(i,:) = ComputeWeights(envelope, thresholds, decRight, snr(i), lenFrame);
end

sigs = categorical({'1.AM', '2.DSB', '3.LSB', '4.USB', '5.LVSB', '6.UVSB', '7.FM', '8.Combined'});
colors = lines(8);

% AM, DSB
figure(8); set(gcf, 'color', 'w');
subplot(1,2,1);
idx = find(decRight == "AM");
hold on;
for i = 1 : sigsNum
    h = bar(sigs(i), weights(idx,i));
    set(h, 'facecolor', colors(i,:));
end
hold off;
grid on; ylim([0, 1]); ylabel('Weight');
title(strcat(decRight(idx), ", SNR = ", num2str(snr(idx)), " dB"));
subplot(1,2,2);
idx = find(decRight == "DSB");
hold on;
for i = 1 : sigsNum
    h = bar(sigs(i), weights(idx,i));
    set(h, 'facecolor', colors(i,:));
end
hold off;
grid on; ylim([0, 1]); ylabel('Weight');
title(strcat(decRight(idx), ", SNR = ", num2str(snr(idx)), " dB"));

% SSB
figure(9); set(gcf, 'color', 'w');
subplot(1,2,1);
idx = find(decRight == "LSB");
hold on;
for i = 1 : sigsNum
    h = bar(sigs(i), weights(idx,i));
    set(h, 'facecolor', colors(i,:));
end
hold off;
grid on; ylim([0, 1]); ylabel('Weight');
title(strcat(decRight(idx), ", SNR = ", num2str(snr(idx)), " dB"));
subplot(1,2,2);
idx = find(decRight == "USB");
hold on;
for i = 1 : sigsNum
    h = bar(sigs(i), weights(idx,i));
    set(h, 'facecolor', colors(i,:));
end
hold off;
grid on; ylim([0, 1]); ylabel('Weight');
title(strcat(decRight(idx), ", SNR = ", num2str(snr(idx)), " dB"));

% VSB
figure(10); set(gcf, 'color', 'w');
subplot(1,2,1);
idx = find(decRight == "LVSB");
hold on;
for i = 1 : sigsNum
    h = bar(sigs(i), weights(idx,i));
    set(h, 'facecolor', colors(i,:));
end
hold off;
grid on; ylim([0, 1]); ylabel('Weight');
title(strcat(decRight(idx), ", SNR = ", num2str(snr(idx)), " dB"));
subplot(1,2,2);
idx = find(decRight == "UVSB");
hold on;
for i = 1 : sigsNum
    h = bar(sigs(i), weights(idx,i));
    set(h, 'facecolor', colors(i,:));
end
hold off;
grid on; ylim([0, 1]); ylabel('Weight');
title(strcat(decRight(idx), ", SNR = ", num2str(snr(idx)), " dB"));
        
% FM, Combined
figure(11); set(gcf, 'color', 'w');
subplot(1,2,1);
idx = find(decRight == "FM");
hold on;
for i = 1 : sigsNum
    h = bar(sigs(i), weights(idx,i));
    set(h, 'facecolor', colors(i,:));
end
hold off;
grid on; ylim([0, 1]); ylabel('Weight');
title(strcat(decRight(idx), ", SNR = ", num2str(snr(idx)), " dB"));
subplot(1,2,2);
idx = find(decRight == "Combined");
hold on;
for i = 1 : sigsNum
    h = bar(sigs(i), weights(idx,i));
    set(h, 'facecolor', colors(i,:));
end
hold off;
grid on; ylim([0, 1]); ylabel('Weight');
title(strcat(decRight(idx), ", SNR = ", num2str(snr(idx)), " dB"));
% % R2017b only
% bar_child = get(bar_h, 'Children');
% mydata=1:sigsNum;
% set(bar_child, 'CData', mydata);
