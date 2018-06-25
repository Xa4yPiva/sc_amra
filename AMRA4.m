function [decision] = AMRA4(keyFeatures, thresholds)
%AMRA4 Summary of this function goes here
%   Detailed explanation goes here

kf = keyFeatures;
thresh = thresholds;
absP = abs(kf.P);

if (absP < thresh.P)
    if (kf.gammaMax < thresh.gammaMax)
        decision = "FM";
    else
        if (kf.sigmaAP < thresh.sigmaAP)
            if (kf.sigmaDP < thresh.sigmaDP)
                decision = "AM";
            else
                decision = "DSB";
            end
        else
            decision = "COmbined";
        end
    end
else
    if (kf.P > 0)
        decision = "LSB";
    else
        if (kf.sigmaAP < thresh.sigmaAP)
            decision = "VSB";
        else
            decision = "USB";
        end
    end
end

end

