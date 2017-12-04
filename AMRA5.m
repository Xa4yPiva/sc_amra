function [decision] = AMRA5(P, yMax, sAP, sDP, threshP, threshY, threshSap, threshSdp)
%AMRA5 Summary of this function goes here
%   Detailed explanation goes here

absP = abs(P);

if (sDP < threshSdp)
    if (absP < threshP)
        decision = "AM";
    else
        decision = "VSB";
    end
else
    if (yMax < threshY)
        decision = "FM";
    else
        if (absP < threshP)
           if (sAP < threshSap)
               decision = "DSB";
           else
               decision = "Combined";
           end
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

