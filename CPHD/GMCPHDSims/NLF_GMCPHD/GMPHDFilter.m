function GMM = GMPHDFilter(GMM, GMM_birth, GMM_spawn, Z, F, Q, H, R, ps, pd, kappa)
%  Function used to implement Gaussian Mixture PHD filter.
%  Inputs: 
%    GMM -- Current Gaussian Mixture Model. This is a structure 
%           with its fields: omega, mean, variance.
%      GMM.omega-- intensity of guassian component.
%      GMM.mean -- mean vector.
%      GMM.variance -- variance matrix.
%    GMM_birth -- Birth target Gaussian Mixture Model. 
%      It's the same tructure as GMM.
%    GMM_spawn -- Spawning target Gaussian Mixture Model. 
%       Same structure as GMM. with extra filed: trans.
%       GMM_spawn.trans -- Transition matrix of spawning target.
%    Z -- Observation set. It's a matrix with each column as a sample.
%    F -- State transition matrix.
%    Q -- Process noise variance matrix.
%    H -- Measurement matrix.
%    R -- Measurement noise variance.
%    ps-- Probability scalar of survival.
%    pd-- Probability scalar of detection.
%    kappa -- Clutter component.
%
%  Outputs:
%    GMM -- The GMM after filtering.
%
%   Ref: Vo B-N, Ma W-K. 
%        The Gaussian mixture probability hypothesis density filter [J]. 
%        Signal Processing, IEEE Transactions on. 2006, 54(11) pp4096 Table I
%
%   Copyright 2012-2013 N.U.D.T. 
%   Author:   David Lieu. 
%   $Revision: 1.0 $  
%   $Date: 2012/03/20 16:42:24 $
%
m = 0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Step1: Prediction for birth and spawning targets.
for n=1:length(GMM_birth)
    m = m + 1;
    WeightPred(m) = GMM_birth(n).omega; %#ok<*AGROW>
    MeanPred(:,m) = GMM_birth(n).mean;
    VarPred(:,:,m) = GMM_birth(n).variance;
end
for n=1:length(GMM_spawn)
    for k=1:length(GMM)
        m = m + 1;
        WeightPred(m) = GMM_spawn(n).omega*GMM(k).omega;
        MeanPred(:,m) = GMM_spawn(n).trans*GMM(k).mean + GMM_spawn(n).mean;
        VarPred(:,:,m) = GMM_spawn(n).trans*GMM(k).variance*...
                        transpose(GMM_spawn(n).trans) +  GMM_spawn(n).variance;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step2: Prediction for existing targets.
for n=1:length(GMM)
    m = m + 1;
    WeightPred(m)  = ps*GMM(n).omega;
    MeanPred(:,m)  = F*GMM(n).mean;
    VarPred(:,:,m) = F*GMM(n).variance*F'+Q;
end
nGMMPred = m;
m = 0;
clear GMM;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step3: Update and output.
for n=1:nGMMPred
    m = m + 1;
    GMM(m).omega    = (1-pd)*WeightPred(n);
    GMM(m).mean     = MeanPred(:,n);
    GMM(m).variance = VarPred(:,:,n);
end
for ell=1:size(Z,2)
    den = kappa(ell);
    idx = m + 1;
    for n=1:nGMMPred
        m     = m + 1;
        zp    = H*MeanPred(:,n);
        S     = H*VarPred(:,:,n)*H' + R;
        Alpha = mvnpdf(Z(:,ell)',zp',S);
        GMM(m).omega = pd*WeightPred(n)*Alpha;
        K = VarPred(:,:,n)*H'/S;
        GMM(m).mean = MeanPred(:,n) + K*(Z(:,ell)- zp);
        GMM(m).variance = (eye(size(K,1))-K*H)*VarPred(:,:,n);
        
        den = den + pd*WeightPred(n)*Alpha;
    end
    for ii=idx:m
        GMM(ii).omega = GMM(ii).omega/den;
    end
end
