function [decision] = AMRA5(keyFeatures, thresholds)
%AMRA5 Summary of this function goes here
%   Detailed explanation goes here

kf = keyFeatures;
thresh = thresholds;
absP = abs(kf.P);

if (kf.sigmaDP < thresh.sigmaDP)
    if (absP < thresh.P)
        decision = "AM";
    else
        decision = "VSB";
    end
else
    if (kf.gammaMax < thresh.gammaMax)
        decision = "FM";
    else
        if (absP < thresh.P)
           if (kf.sigmaAP < thresh.sigmaAP)
               decision = "DSB";
           else
               decision = "Combined";
           end
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

