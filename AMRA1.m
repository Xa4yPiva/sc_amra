function [decision] = AMRA1(keyFeatures, thresholds)
%AMRA1 Summary of this function goes here
%   Detailed explanation goes here

kf = keyFeatures;
thresh = thresholds;
absP = abs(kf.P);

if (kf.sigmaDP < thresh.sigmaDP)
    if (absP < thresh.P / 2)
        decision = "AM";
    else
        if (kf.P > 0)
            decision = "LVSB";
        else
            decision = "UVSB";
        end
    end
else
    if (absP < thresh.P)
        if (kf.gammaMax < thresh.gammaMax)
            if (kf.sigmaAP > thresh.sigmaAP)
                decision = "FM";
            else
                decision = "Unknown";
            end
        else
            if (kf.sigmaAP < thresh.sigmaAP)
                decision = "DSB";
            else
                decision = "Combined";
            end
        end
    else
        if (kf.P > 0)
            decision = "LSB";
        else
            decision = "USB";
        end
    end
end

end

