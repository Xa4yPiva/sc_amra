function [ weights ] = ComputeWeights( envelope, thresholds, decRight, SNR, lenFrame )
%COMPUTEWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

sigsNum = length(decRight);
lenSignal = length(envelope);
framesNum = floor(lenSignal / lenFrame);
w = zeros(1, sigsNum);
cyclesNum = framesNum;
iteration = 0;
h = waitbar(0, 'Computing weights...');
for i = 1 : framesNum
    env = envelope((i-1)*lenFrame+1 : i*lenFrame);
    env = awgn(env, SNR, 'measured');
    kf = KeyFeatures(env, thresholds.ampl);
    decision = AMRA1(kf, thresholds);
    idx = find(decRight == decision);
    w(idx) = w(idx) + 1;
    iteration = iteration + 1;
    waitbar(iteration / cyclesNum);
end
close(h);
weights = w / framesNum;

end

