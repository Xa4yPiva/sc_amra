function [decision] = AMRA4(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp)
%AMRA4 Summary of this function goes here
%   Detailed explanation goes here

absP = abs(P);

if (absP < threshP)
    if (yMax < threshY)
        decision = "FM";
    else
        if (sAP < threshSap)
            if (sDP < threshSdp)
                decision = "AM";
            else
                decision = "DSB";
            end
        else
            decision = "COmbined";
        end
    end
else
    if (P > 0)
        decision = "LSB";
    else
        if (sAP < threshSap)
            decision = "VSB";
        else
            decision = "USB";
        end
    end
end

end

