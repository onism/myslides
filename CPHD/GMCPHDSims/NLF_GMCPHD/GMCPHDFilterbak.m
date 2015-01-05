function [GMM,Prob]= GMCPHDFilter(GMM, Prob, GMM_birth, ProbBirth, Z, F, Q, H, R, ps, pd, Clutters)
%  Function used to implement Gaussian Mixture Cardinality PHD filter.
%  Inputs: 
%    GMM -- Current Gaussian Mixture Model. This is a structure 
%           with its fields: omega, mean, variance.
%      GMM.omega-- intensity of guassian component.
%      GMM.mean -- mean vector.
%      GMM.variance -- variance matrix.
%    GMM_birth -- Birth target Gaussian Mixture Model. 
%      It's the same tructure as GMM.
%    Prob -- Cardinality distribution.
%    ProbBirth -- Cardinality distribution of birth target.
%    Z -- Observation set. It's a matrix with each column as a sample.
%    F -- State transition matrix.
%    Q -- Process noise variance matrix.
%    H -- Measurement matrix.
%    R -- Measurement noise variance.
%    ps-- Probability scalar of survival.
%    pd-- Probability scalar of detection.
%    Clutters -- Clutter cardinality distribution and spatial distribution.
%      Clutters.funCardPdf -- function handle used to specify the cardinality distribution.
%      Clutters.funSpatialDist -- function handle used to specify the spatial distribution of
%      clutter.
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
global JJgl;
global KKgl;
global NNgl;
global WWgl;
global LLgl;

m     = 0;
nmax  = length(Prob)-1;
kappa = Clutters.funCardPdf;
c     = Clutters.funSpatialDist;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Step1: Prediction for birth and spawning targets.
for n=1:length(GMM_birth)
    m = m + 1;
    WeightPred(m) = GMM_birth(n).omega; %#ok<*AGROW>
    MeanPred(:,m) = GMM_birth(n).mean;
    VarPred(:,:,m) = GMM_birth(n).variance;
end
% for n=1:length(GMM_spawn)
%     for k=1:length(GMM)
%         m = m + 1;
%         WeightPred(m) = GMM_spawn(n).omega*GMM(k).omega;
%         MeanPred(:,m) = GMM_spawn(n).trans*GMM(k).mean + GMM_spawn(n).mean;
%         VarPred(:,:,m) = GMM_spawn(n).trans*GMM(k).variance*...
%                         transpose(GMM_spawn(n).trans) +  GMM_spawn(n).variance;
%     end
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step2: Prediction for existing targets.
for n=1:length(GMM)
    m = m + 1;
    WeightPred(m) = ps*GMM(n).omega;
    MeanPred(:,m) = F*GMM(n).mean;
    VarPred(:,:,m) = F*GMM(n).variance*F'+Q;
end
nGMMPred = m;
m        = 0;
clear GMM;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step3: Prediction of Cardinality distribution.
ProbPred = sum(WWgl.*(1-ps).^JJgl.*...
               Prob(JJgl+1).*...
               (ps/(1-ps)).^KKgl.*...
               ProbBirth(NNgl-KKgl+1).*...
               double(LLgl>0),1);

% % % ProbPred = zeros(size(Prob));
% % % UIdx     = triu(repmat(1:nmax+1,[nmax+1,1]),0);
% % % RIdx     = triu(repmat((1:nmax+1)',[1,nmax+1]),0);

% % for n=0:nmax   
% % %     MatIdx = UIdx(1:n+1,:)';  % Use matrix operation for speed.
% % %     JJ     = MatIdx(MatIdx>0);
% % %     JJ     = JJ' - 1;
% % %     
% % %     MatIdx = RIdx(1:n+1,:)';
% % %     KK     = MatIdx(MatIdx>0);
% % %     KK     = KK' - 1;
% % %     ProbPred(n+1) = sum(factorial(JJ)./...
% % %                         factorial(JJ-KK)./...
% % %                         factorial(KK).*...
% % %                         (1-ps).^(JJ-KK).*...
% % %                         Prob(JJ+1).*...
% % %                         ps.^KK.*ProbBirth(n-KK+1));
%  For Loop with slow mathod               
% %     for kk=0:n  
% %         hbar = sum(factorial(kk:nmax)./factorial((kk:nmax)-kk)/factorial(kk).*...
% %                (1-ps).^((kk:nmax)-kk).*Prob((kk:nmax)+1));
% %         ProbPred(n+1) = ProbPred(n+1) + ps^kk*ProbBirth(n-kk+1)*hbar;
% %     end
% % end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step4: Update the Cardinality Distribution.
nz         = size(Z,2);
Likelihood = zeros(nGMMPred,nz);

for ii=1:nGMMPred
    zp   = H*MeanPred(:,ii);
    S    = H*VarPred(:,:,ii)*H' + R;
    for jj=1:nz
        Likelihood(ii,jj) = mvnpdf(Z(:,jj)',zp',S)/c(Z(:,jj));
    end
end

%% Note that: This is not the direct way to calculate the weight series.
%   Matrix operations are used to avoid "for loops" in order to implement quick calculations.
%  The direct routine follow the reference is implemented in the following lines
%  in comments. 
if nmax>=nz
    UIdx = triu(repmat((1:nz+1)',[1,nz+1]),0);
    UIdx = [UIdx repmat((1:nz+1)',[1,nmax-nz])];
    NN   = repmat(1:nmax+1,[nz+1,1])-1;    
else
    UIdx = triu(repmat((1:nmax+1)',[1,nmax+1]),0);
    NN   = repmat(1:nmax+1,[nmax+1,1])-1;
end


JJ = UIdx - 1;
JJ(JJ==-1) = 0;  % in case the factorial error with negative inputs.
KK = NN - JJ;
KK(KK<0) = 0;    % in case the factorial error with negative inputs.
KK1 = NN - JJ -1;
KK1(KK1<0) = 0;

WW = factorial(nz-JJ).*(pd/(1-pd)).^JJ.*(1-pd).^NN.*...
      ElemSymmFunc((WeightPred*Likelihood)',JJ).*...
      kappa(nz-JJ)./sum(WeightPred).^JJ.*factorial(NN);
Weight0 = sum(WW./factorial(KK).*double(UIdx>0),1);
Weight1 = sum(WW./factorial(KK1)/...
              (1-pd)/sum(WeightPred).*triu(double(UIdx>0),1),1);

%% Direct implementation with "for loops".
% for n=0:nmax
%%  This is the optimized code for fast calculation.
%     jj=0:min([nz,n]);
%     WW = factorial(nz-jj).*pd.^jj.*(1-pd).^(n-jj).*...
%          ElemSymmFunc((WeightPred*Likelihood)',jj).*...
%          kappa(nz-jj)./sum(WeightPred).^jj;
%     Weight0(n+1) = sum(WW.*factorial(jj).*ElemSymmFunc(ones(n,1), jj));
%     Weight1(n+1) = sum(WW.*factorial(jj+1).*ElemSymmFunc(ones(n,1), jj+1)/...
%          (1-pd)/sum(WeightPred));

%%%%% This is the orginal code for the reader to undertand the CPHD weight calculation.
%     Weight0(n+1) = sum(factorial(nz-jj).*factorial(jj).*......
%         ElemSymmFunc(ones(n,1), jj).*...
%         pd.^jj.*(1-pd).^(n-jj).*...
%         ElemSymmFunc((WeightPred*Likelihood)',jj).*...
%         kappa(nz-jj)./...
%         sum(WeightPred).^jj);
%     Weight1(n+1) = sum(factorial(nz-jj).*factorial(jj+1).*...
%         ElemSymmFunc(ones(n,1), jj+1).*...
%         pd.^jj.*(1-pd).^(n-jj-1).*...
%         ElemSymmFunc((WeightPred*Likelihood)',jj).*...
%         kappa(nz-jj)./...
%         sum(WeightPred).^(jj+1));
    
%     for ii=1:nz        
%         LMatrix = [Likelihood(:,1:ii-1) Likelihood(:,ii+1:end)];
%         jj=0:min([nz-1,n]);
%         WeightZ(n+1, ii) =sum(factorial(nz-1-jj).*factorial(jj+1).*...
%               ElemSymmFunc(ones(n,1), jj+1).*...
%               pd.^jj.*(1-pd).^(n-jj-1).*...
%               ElemSymmFunc((WeightPred*LMatrix)',jj).*...
%               kappa(nz-1-jj)./...
%               sum(WeightPred).^(jj+1));
%     end
% end
Prob = (Weight0.*ProbPred);
Prob = Prob/sum(Prob);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step5: Update PHD and output.
for n=1:nGMMPred
    m = m + 1;
    GMM(m).omega    = sum(Weight1.*ProbPred)/sum(Weight0.*ProbPred)*...
                       (1-pd)*WeightPred(n);
    GMM(m).mean     = MeanPred(:,n);
    GMM(m).variance = VarPred(:,:,n);
end


% Calculating the weight matrix exclude current observation.
% Note that: This is not the direct way to calculate the weight series.
%     I use the matrix operations to implement quick calculations.
if nmax>=nz-1
    UIdx = triu(repmat((1:nz)',[1,nz]),0);
    UIdx = [UIdx repmat((1:nz)',[1,nmax-nz+1])];
    NN   = repmat(1:nmax+1,[nz,1])-1;    
else
    UIdx = triu(repmat((1:nmax+1)',[1,nmax+1]),0);
    NN   = repmat(1:nmax+1,[nmax+1,1])-1;
end
JJ = UIdx - 1;
JJ(JJ==-1) = 0;  % in case the factorial function error with negative inputs.
KK1 = NN - JJ - 1;
KK1(KK1<0) = 0;  % in case the factorial function error with negative inputs.
WW0  = factorial(nz-1-JJ)./factorial(KK1).*...% .*factorial(NN)
       pd.^JJ.*(1-pd).^KK1.*kappa(nz-1-JJ)./...                  
       sum(WeightPred).^(JJ+1).*triu(double(UIdx>0),1);
for ell=1:size(Z,2)
    LMatrix = [Likelihood(:,1:ell-1) Likelihood(:,ell+1:end)];
    WW      = sum(ElemSymmFunc((WeightPred*LMatrix)',JJ).*WW0,1);
    for n=1:nGMMPred
        m = m + 1;
        zp = H*MeanPred(:,n);
        S = H*VarPred(:,:,n)*H' + R;
%         Alpha = mvnpdf(Z(:,ell)',zp',S)*...
%             sum(WeightZ(:,ell)'.*ProbPred)/sum(Weight0.*ProbPred)/c(Z(:,ell));
        Alpha = mvnpdf(Z(:,ell)',zp',S)*...
                sum(WW.*ProbPred)/sum(Weight0.*ProbPred)/c(Z(:,ell));
        
        K = VarPred(:,:,n)*H'/S;
        GMM(m).omega    = pd*WeightPred(n)*Alpha;
        GMM(m).mean     = MeanPred(:,n) + K*(Z(:,ell)- zp);
        GMM(m).variance = (eye(size(K,1))-K*H)*VarPred(:,:,n);
    end
end
