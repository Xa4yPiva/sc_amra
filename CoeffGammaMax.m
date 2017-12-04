function gammaMax = CoeffGammaMax(aCentNorm, Nfft)
%COEFFGAMMAMAX Summary of this function goes here
%   Detailed explanation goes here

aAbsFft = fftshift(abs(fft(aCentNorm, Nfft))); % Check this ( /N or /(N/2) or nothing)
gammaMax = max(aAbsFft .* aAbsFft / Nfft);

end

