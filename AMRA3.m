function [decision] = AMRA3(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp)
%AMRA3 Summary of this function goes here
%   Detailed explanation goes here

absP = abs(P);

if (yMax < threshY)
    decision = "FM";
else
    if (sAP < threshSap)
        if (absP < threshP)
            if (sDP < threshSdp)
                decision = "AM";
            else
                decision = "DSB";
            end
        else
            decision = "VSB";
        end
    else
        if (absP < threshP)
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

