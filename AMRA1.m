function [decision] = AMRA1(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp)
%AMRA1 Summary of this function goes here
%   Detailed explanation goes here

absP = abs(P);

if (sDP < threshSdp)
    if (absP < threshP)
        decision = "AM";
    else
        decision = "VSB";
    end
else
    if (absP < threshP)
        if (yMax < threshY)
            decision = "FM";
        else
            if (sAP < threshSap)
                decision = "DSB";
            else
                decision = "Combined";
            end
        end
    else
        if (P > 0)
            decision = "LSB";
        else
            decision = "USB";
        end
    end
end

end

