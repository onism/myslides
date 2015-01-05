function GMM = GMPHDPrunning(GMM, Jmax, T, U)
% Routine used to truncate these gaussian componemts with small weight.
%  And merge these components with almost the same means.
%
% Inputs:
%   GMM -- Struct which is going to be truncate and unite.
%   Jmax -- maximum allowable number of Gaussian terms.
%   T  -- truncation threshold.
%   U  -- merging threshold.
%
% Outputs:
%   GMM -- guassian mixture model after truncation and merging.
%
%   Ref: Vo B-N, Ma W-K. 
%        The Gaussian mixture probability hypothesis density filter [J]. 
%        Signal Processing, IEEE Transactions on. 2006, 54(11) pp4097 Table II
%
%   Copyright 2012-2013 N.U.D.T. 
%   Author:   David Lieu. 
%   $Revision: 1.0 $  
%   $Date: 2012/03/20 16:42:24 $
%
%
cnt = 0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step1: truncate these Gaussian Component with small weight.
if isempty(GMM)
    return;
end
for n=1:length(GMM)
    if GMM(n).omega > T
        cnt = cnt + 1;
        II(cnt) = n; %#ok<*AGROW,*NASGU>
    end    
end
cnt = 0;
nDimX = size(GMM(1).mean,1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step2: merge these Gaussian Components with same means.
while ~isempty(II)
    cnt     = cnt + 1;
    Weights = [GMM(II).omega];
    [~,idx] = max(Weights);    
    w = 0;
    X = zeros(nDimX,1);
    P = zeros(nDimX);
    w0 = GMM(II(idx)).omega;
    X0 = GMM(II(idx)).mean;
    P0 = GMM(II(idx)).variance;
    
    for n=1:length(Weights)
        d = (GMM(II(n)).mean-X0)'/GMM(II(n)).variance*(GMM(II(n)).mean-X0); %#ok<*MINV>
        if d<=U
            w = w + GMM(II(n)).omega;
            X = X + GMM(II(n)).omega*GMM(II(n)).mean;
            P = P + GMM(II(n)).omega*(GMM(II(n)).variance + ...
                   (GMM(II(n)).mean-X0)*(GMM(II(n)).mean-X0)');
            II(n) = 0;
        end
    end
%     X = X/w;
%     P = P/w;
    GMMPruning(cnt).omega = w;
    GMMPruning(cnt).mean = X/w;
    GMMPruning(cnt).variance = P/w;
    
    II = II(II>0);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Step3: If the number of  Gaussian components exceed the maximum allowable numbers
%     Choose those who have the largest weights.
if cnt>Jmax
    Weights = [GMMPruning.omega];
    [~, idx] = sort(Weights, 'descend');
    GMM = GMMPruning(idx);
else
    GMM = GMMPruning;
end