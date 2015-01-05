function Xhat = GMPHDStateExtractor(GMM)
% Routine used to get target state from Gaussian Mixture Components
%
% Inputs:
%   GMM -- Struct which is going to be truncate and unite.
%
% Outputs:
%   Xhat -- Multiple target state.
%
%   Ref: Vo B-N, Ma W-K. 
%        The Gaussian mixture probability hypothesis density filter [J]. 
%        Signal Processing, IEEE Transactions on. 2006, 54(11) pp4097 Table III
%
%   Copyright 2012-2013 N.U.D.T. 
%   Author:   David Lieu. 
%   $Revision: 1.0 $  
%   $Date: 2012/03/20 16:42:24 $
%
%
Weights = [GMM.omega];
Xhat = [];
for ii=1:length(Weights)
    if Weights(ii)>0.5
        for jj=1:round(Weights(ii))
            Xhat = [Xhat GMM(ii).mean]; %#ok<*AGROW>
        end
    end
end