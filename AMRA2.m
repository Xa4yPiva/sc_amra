function [decision] = AMRA2(keyFeatures, thresholds)
%AMRA2 Summary of this function goes here
%   Detailed explanation goes here

kf = keyFeatures;
thresh = thresholds;
absP = abs(kf.P);

if (kf.gammaMax < thresh.gammaMax)
   decision = "FM";
else
    if (kf.sigmaAP < thresh.sigmaAP)
        if (kf.sigmaDP < thresh.sigmaDP)
            if (absP < thresh.P)
                decision = "AM";
            else
                decision = "VSB";
            end
        else
            decision = "DSB";
        end
    else
        if (absP < thresh.P)
            decision = "Combined";
        else
            if (P < 0)
                decision = "USB";
            else
                decision = "LSB";
            end
        end
    end
end
    
end

