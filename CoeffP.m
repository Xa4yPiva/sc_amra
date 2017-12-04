function P = CoeffP(envelope, Nfft)
%COEFFP Summary of this function goes here
%   Detailed explanation goes here

spectrumAbs = abs(fftshift(fft(envelope, Nfft)));
fcn = Nfft / 2 + 1; % Check this
Pl = sum(spectrumAbs(1 : fcn)   .* spectrumAbs(1 : fcn));
Pu = sum(spectrumAbs(fcn : end) .* spectrumAbs(fcn : end));
P = (Pl - Pu) / (Pl + Pu);

end

