function [decision] = AMRA_Khm(keyFeatures, thresholds)
%AMRA_KHM Summary of this function goes here
%   Detailed explanation goes here

kf = keyFeatures;
thresh = thresholds;

% if (kf.sigmaDP < thresh.sigmaDP)
%     decision = "AM";
% else
%     if (absP < thresh.P)
%         if (kf.gammaMax < thresh.gammaMax)
%             if (kf.sigmaAP > thresh.sigmaAP)
%                 decision = "FM";
%             else
%                 decision = "Unknown";
%             end
%         else
%             if (kf.sigmaAP < thresh.sigmaAP)
%                 decision = "DSB";
%             else
%                 decision = "Combined";
%             end
%         end
%     else
%         if (kf.P > 0)
%             decision = "LSB";
%         else
%             decision = "USB";
%         end
%     end
% end

end

