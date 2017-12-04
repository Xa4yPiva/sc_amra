function [decision] = AMRA2(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp)
%AMRA2 Summary of this function goes here
%   Detailed explanation goes here

absP = abs(P);

if (yMax < threshY)
   decision = "FM";
else
    if (sAP < threshSap)
        if (sDP < threshSdp)
            if (absP < threshP)
                decision = "AM";
            else
                decision = "VSB";
            end
        else
            decision = "DSB";
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

