function P = CoeffP(signal, Nfft, fc, fs)
%COEFFP Summary of this function goes here
%   Detailed explanation goes here

specSignal = fftshift(fft(signal, Nfft));
specOneSide = specSignal(end/2+1 : end);
Nfft = length(specOneSide) * 2;
fcn = fc * Nfft / fs + 1; % Check this
Pl = sum(abs(specOneSide(1 : fcn))   .* abs(specOneSide(1 : fcn)));
Pu = sum(abs(specOneSide(fcn : end)) .* abs(specOneSide(fcn : end)));
P = (Pl - Pu) / (Pl + Pu);

end

