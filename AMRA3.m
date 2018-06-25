function [decision] = AMRA3(keyFeatures, thresholds)
%AMRA3 Summary of this function goes here
%   Detailed explanation goes here

kf = keyFeatures;
thresh = thresholds;
absP = abs(kf.P);

if (kf.gammaMax < thresh.gammaMax)
    decision = "FM";
else
    if (kf.sigmaAP < thresh.sigmaAP)
        if (absP < thresh.P)
            if (kf.sigmaDP < thresh.sigmaDP)
                decision = "AM";
            else
                decision = "DSB";
            end
        else
            decision = "VSB";
        end
    else
        if (absP < thresh.P)
            decision = "Combined";
        else
            if (kf.P < 0)
                decision = "USB";
            else
                decision = "LSB";
            end
        end
    end
end

end

