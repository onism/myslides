% Scripts used to implement the Gaussian Mixture Filter.
%  Indivisual target motion model and observation model both are linear and gaussian.
%
%
%   Ref: Vo B-N, Ma W-K. 
%        Analytic Implementations of the Cardinalized Probability Hypothesis Density Filter[J]. 
%        Signal Processing, IEEE Transactions on. 2007, 55(7) pp2559 Section IV
%
%   Remark: This is not the MAIN File.
%
%   Copyright 2012-2013 N.U.D.T. 
%   Author:   David Lieu. 
%   $Revision: 1.0 $  
%   $Date: 2012/03/20 16:42:24 $
%
%
if ~exist('MATLAB_DEBUG','var')
    MATLAB_DEBUG = 1;
end
if MATLAB_DEBUG
    close all; clear;
    SystemConstDeclare;
    GMCPHDInitialParameters;
end

%% Variables preallocations for speed.
ZCLobsv   = cell(nSims, nMCs);  % Observation cells.
XCardTrue = zeros(nSims, nMCs); % True cardinality of multiple targets.
XCardHat  = zeros(nSims, nMCs); % Target cardinality estimation.
SeqGMM    = cell(nSims, nMCs);  % Sequential Gaussian Mixture Model.
XCLFilter = cell(nSims, nMCs);
OSPA      = zeros(nSims, nMCs); % Optimal SubPattern Assignment Metric.
Hausdorf  = OSPA;   % Hausdorf Metric between true RFS and Filtered RFS.
OMAT      = zeros(nSims, nMCs, 2);   % Optimal MAss Assignment Transfer Metric.

GMM(1).omega = 0;
GMM(1).mean  = zeros(nDimX,1);
GMM(1).variance = zeros(nDimX);

%% Set the clutter cardinality and intensity distribution.
Clutters.funCardPdf = @(x) poisspdf(x,lambdac*vol); % here cardinality distribution is poisson.
Clutters.funSpatialDist = @(x) 1/vol; %  the spatial distribution is uniform.

%% Initialize the progress bar.
% if ~MATLAB_DEBUG || nMCs>1,
%     flagShowWaitbar = 1;
%     hWaitbar = waitbar(0, ['Monte Carlo Running... lambdac=', num2str(lambdac)]);
% end
%% Monte Carlo runs...
for m=1:nMCs
%     if flagShowWaitbar
%         waitbar(m/nMCs);
%     end
    k = 1;
    % SeqGMM{k,m}    = [];
    SeqGMM{k,m}    = GMM_birth;
    CardDist{k,m}  = [ones(1,4)/4 zeros(1, nCardMax-4)];
    XCLFilter{k,m} = [GMM_birth.mean];
    XCardHat(k,m)  = 0;
    
    X = [];
    for n=1:length(Track)
        if Track(n).isactive(m)
            Track(n).len(m) = Track(n).len(m) + 1;
            if ~flagSpecified
                offset = Track(n).offset;
                switch Track(n).type
                    case TYPE_BIRTH
                        Track(n).state(:,k,m) = mvnrnd(GMM_birth(offset).mean,...
                             GMM_birth(offset).variance)';                        
                    case TYPE_SPAWN
                        Track(n).state(:,k,m) = mvnrnd(GMM_spawn(offset).mean,...
                             GMM_spawn(offset).variance)';
                end
                
            end
            X = [X, Track(n).state(:,k,m)];
        end
    end
        
    XCardTrue(1,m) = 3;
    
    Z = [];
    for n=1:size(X,2)
        if rand < probdt,  % for probability of detection
            Z = [Z, mvnrnd(H*(X(:,n)), R)']; %#ok<*AGROW>
        end
    end
    nc = poissrnd(lambdac*vol);    
    ZCLobsv{1,m} = [Z [xmin + (xmax - xmin)*rand(1,nc);...
                       ymin + (ymax - ymin)*rand(1,nc)]];
    
    for k=1:nSims-1
        %% For existing targets.
%         if ~flagShowWaitbar
%             k %#ok<*NOPTS>
%         end
        disp(['>>> Running in m=' num2str(m) '/' num2str(nMCs) ', k=' num2str(k) '/' num2str(nSims) '...']);
        
        Xtrue = [];
        for n=1:length(Track)
            Track(n).isactive(m) = Track(n).isactive(m) &&...
                     (k>=Track(n).t0(m)) && (k<=Track(n).tf(m));
            if Track(n).isactive(m)  % Current track is active?
                if rand < probsv
                    X = Track(n).state(:,k,m);
                                           
                    switch Track(n).model(k)
                        case MODEL_CV  % Just update the possition part of State vector.
                            X(1:4) = Fcv*X(1:4);
                            %X(1:4) = mvnrnd((Fcv*X(1:4))',Qcv)';
                        case MODEL_CT  % Update full state vector.
                            X = Fct(X(5),T,alpha_man)*X;
                            %X = mvnrnd((Fct(X(5),T,alpha_man)*X)', Qct)';
                    end
                
                    
                    flagInFOV = (X(1)>=xmin) && (X(1)<=xmax);
                    flagInFOV = flagInFOV && (X(3)>=ymin) && (X(3)<=ymax);
                                        
                    if flagInFOV  % Save the true track.
                        Track(n).state(1:4,k+1,m) = X(1:4);
                        Track(n).len(m) = Track(n).len(m) + 1; %#ok<*SAGROW>
                        Xtrue = [Xtrue X];
                        XCardTrue(k+1,m) = XCardTrue(k+1,m) + 1;
                    else
                        Track(n).tf(m) = k;
                        Track(n).isactive(m) = 0;                        
                    end
                else % Current track died.
                    Track(n).tf(m) = k;
                    Track(n).isactive(m) = 0;
                end
            else % current track is not alive.
                if ~Track(n).len(m) % Current track has never been updated before. So it must be track to be activated.
                    if k==Track(n).t0(m)-1 % activate current track at next time step.
                        Track(n).isactive(m) = 1;                        
                    end
                end
            end
        end
       
        %% Generate the observation set.
        Z = [];
        for n=1:size(Xtrue,2)
            if rand < probdt,  % for probability of detection
                Z = [Z, mvnrnd((H*(Xtrue(:,n)))', R)'];
            end
        end
        nc = poissrnd(lambdac*vol);
        ZCLobsv{k+1,m} = [Z [xmin + (xmax-xmin)*rand(1,nc);...
                             ymin + (ymax-ymin)*rand(1,nc)]];
        
        %% Start PHD Filtering Routine.
        Z = ZCLobsv{k+1,m}; 
%        ProbBirth = poisspdf(0:nCardMax,8);
%        ProbBirth = [ones(1,8)/8 zeros(1, nCardMax-8)];
        SeqGMM{k+1,m} = GMPHDFilter(SeqGMM{k,m}, GMM_birth, GMM_spawn, Z, F, Q, H, R, pspar, pdpar,...
                        ones(size(Z,2),1)*lambdac*vol*u);
%         [SeqGMM{k+1,m}, CardDist{k+1,m}] = GMCPHDFilter(SeqGMM{k,m}, CardDist{k,m}, GMM_birth,...
%             ProbBirth, Z, F, Q, H, R, pspar, pdpar, Clutters);  
        SeqGMM{k+1,m} = GMPHDPrunning(SeqGMM{k+1,m}, nGMMax, thtrunc, thunion);
        %% Multiple Target State Extraction.
        XCLFilter{k+1,m} = GMPHDStateExtractor(SeqGMM{k+1,m});
        
        %% Target RFS cardinality estimation.
        XCardPred = pspar;
        for n=1:length(GMM_spawn)
            XCardPred = XCardPred + GMM_spawn(n).omega;
        end
        XCardPred = XCardHat(k,m)*XCardPred;
        for n=1:length(GMM_birth)
            XCardPred = XCardPred + GMM_birth(n).omega;
        end
        XCardHat(k+1, m) = (1-pdpar)*XCardPred;
        for n=1:length(SeqGMM{k+1,m})
            GMM = SeqGMM{k+1, m};
            XCardHat(k+1, m) = XCardHat(k+1, m) + GMM(n).omega;
        end
%         [~,XCardHat(k+1,m)] = max(CardDist{k+1,m});
%         XCardHat(k+1,m)     = XCardHat(k+1,m) - 1;
        
        %% Calculate the OSPA and Metrics between True State RFS and Filtered RFS.
        Xhat = [];
        Xt   = [];
        if ~isempty(XCLFilter{k+1,m})
            Xhat = XCLFilter{k+1,m};%([1,3],:);
        end
        if ~isempty(Xtrue)
            Xt = Xtrue;%([1 3],:);
        end
        OSPA(k+1, m) = ospa_dist(Xhat,Xt,100,2);
        Hausdorf(k+1,m) = Hausdorf_dist(Xhat,Xt);
        OMAT(k+1,m,1) = omat_dist(Xhat, Xt, 2);
        OMAT(k+1,m,2) = omat_dist(Xhat, Xt, 3);
    end    
end
% if flagShowWaitbar
%     close(hWaitbar);
% end
Hlg = [1 0 0 0 0; 0 0 1 0 0];
[CPEP, CardPosExpect] = CalcTrackMeasure(Track, XCLFilter, H, cpep_radius);
% CPEP -- Circular position error probability.
% CardPosExpect -- Expected absolute error on target number.
%% Figure results.
m = 1;
nMarker = length(STR_MARKER);
figure(1); 
subplot(211);whitebg('white'); cla; hold on;
for k=1:nSims
    Z = ZCLobsv{k,m};
    if ~isempty(Z)
        plot(k*ones(1,size(Z,2)), Z(1,:), 'kx',...
            'MarkerSize',6,'color', 0.85*[1 1 1]);
    end
    X = XCLFilter{k,m};
    if ~isempty(X)
        plot(k*ones(1,size(X,2)), X(1,:), 'o', 'color',...
                      0.3*[1 1 1],'MarkerSize',4);
    end
end
for n=1:length(Track)
    if Track(n).len
        t = Track(n).t0(m):Track(n).tf(m);
        Xtrue = Track(n).state(:,t,m);
        hTargetFigure(n) = ...
            plot(t, Xtrue(1,:),'k-','linewidth',2);
%             plot(t, Xtrue(1,:),['k-' STR_MARKER(mod(n,nMarker)+1)],...
%                                 'linewidth',2.5,'MarkerSize',4);        
        strLegend{n} = Track(n).id;
    end
end
hold off;
xlabel('time step');
ylabel('X Position(m)');
title('Trajactory of Multiple Objects in X axis.');
% hLegend = legend(hTargetFigure, strLegend,...
%             'location','Best', 'Orientation','horizontal');
% set(hLegend, 'box', 'on');
set(gca, 'Box', 'on');
axis([1 nSims xmin xmax]);

subplot(212);whitebg('white'); cla; hold on;
for k=1:nSims
    Z = ZCLobsv{k,m};
    if ~isempty(Z)
        plot(k*ones(1,size(Z,2)), Z(2,:), 'kx',...
            'MarkerSize',6,'color', 0.85*[1 1 1]);
    end
    X = XCLFilter{k,m};
    if ~isempty(X)
        plot(k*ones(1,size(X,2)), X(3,:), 'o', 'color',...
                      0.3*[1 1 1],'MarkerSize',4);
    end
end
for n=1:length(Track)
    if Track(n).len
        t = Track(n).t0(m):Track(n).tf(m);
        Xtrue = Track(n).state(:,t,m);   
          plot(t, Xtrue(3,:),'k-','linewidth',2);
%         plot(t, Xtrue(3,:),['k-' STR_MARKER(mod(n,nMarker)+1)],...
%                                 'linewidth',2.5,'MarkerSize',4);
    end
end
hold off
axis([1 nSims ymin ymax]);
xlabel('time step');
ylabel('Y Position(m)');
title('Trajactory of Multiple Objects in Y axis.');
set(gca, 'Box', 'on');
% print -depsc GMM_Trajactory.eps
% saveas(gcf, 'Trajactory', 'fig');

figure(2);
whitebg('white'); cla; hold on;
for n=1:length(Track)
    if Track(n).len
        t = Track(n).t0(m):Track(n).tf(m);
        Xtrue = Track(n).state(:,t,m);
        plot(Xtrue(1,:), Xtrue(3,:), 'k-');
        hAx(1) = plot(Xtrue(1,1), Xtrue(3,1), 'ko');
        hAx(2) = plot(Xtrue(1,end), Xtrue(3,end), 'ks');
    end
end
for k=1:nSims
    X = XCLFilter{k,m};
    if ~isempty(X)
        plot(X(1,:), X(3,:), 'o', 'color',...
                      0.3*[1 1 1],'MarkerSize',4);
    end
end
xlabel('X Position(m)');
ylabel('Y Position(m)');
title('Trajactory of Multiple Objects in Plane.');
axis equal;
set(gca, 'Box', 'on');
legend(hAx, 'Start Position', 'End Position', 'location', 'Best');
% print -depsc GMM_TrajactoryInPlane.eps
hold off
axis([xmin xmax ymin ymax]);

figure(3)
whitebg('white'); cla; hold on;
plot(1:nSims, XCardTrue(:,m), 'k-','linewidth',2.5)
plot(2:nSims, round(XCardHat(2:end, m)), 'ko');
xlabel('time step');
ylabel('n_T');
title('True cardinality VS Estimation');
legend('True', 'Estimation');
set(gca, 'Box', 'on');
hold off;

figure(4)
subplot(211);whitebg('white'); cla; hold on;
for k=1:nSims
    X = XCLFilter{k,m};
    if ~isempty(X)
        plot(k*ones(1,size(X,2)), X(2,:), 'o', 'color',...
                      0.3*[1 1 1],'MarkerSize',4);
    end
end
for n=1:length(Track)
    if Track(n).len
        t = Track(n).t0(m):Track(n).tf(m);
        Xtrue = Track(n).state(:,t,m);
        hTargetFigure(n) = ...
            plot(t, Xtrue(2,:),['k-' STR_MARKER(mod(n,nMarker)+1)],...
                                'linewidth',2.5,'MarkerSize',4);
        %strLegend{n} = Track(n).id;
    end
end

xlabel('time step');
ylabel('X Velocity(m/s)');
title('Velocity of Multiple Objects in X axis.');
% hLegend = legend(hTargetFigure, strLegend,...
%             'location','Best', 'Orientation','horizontal');
% set(hLegend, 'box', 'on');
set(gca, 'Box', 'on');
hold off;

subplot(212);whitebg('white'); cla; hold on;
for k=1:nSims
    X = XCLFilter{k,m};
    if ~isempty(X)
        plot(k*ones(1,size(X,2)), X(4,:), 'o', 'color',...
                      0.3*[1 1 1],'MarkerSize',4);
    end
end
for n=1:length(Track)
    if Track(n).len
        t = Track(n).t0(m):Track(n).tf(m);
        Xtrue = Track(n).state(:,t,m);        
        plot(t, Xtrue(4,:),['k-' STR_MARKER(mod(n,nMarker)+1)],...
                                'linewidth',2.5,'MarkerSize',4);
    end
end
xlabel('time step');
ylabel('Y Velocity(m/s)');
title('Trajactory of Multiple Objects in Y axis.');
set(gca, 'Box', 'on');
hold off;

figure(5);
subplot(211);whitebg('white'); cla; 
plot(2:nSims, CPEP(2:end), 'k-','linewidth',2.5);
hold on;
plot(0:nSims, mean(CPEP)*ones(1,nSims+1), 'k-');
xlabel('time step');
ylabel('CPEP');
title('Circular Position Error Probability VS. Time step');
set(gca, 'Box', 'on');
subplot(212);whitebg('white'); cla
plot(2:nSims, CardPosExpect(2:end), 'k-','linewidth',2.5);
hold on;
plot(0:nSims, mean(CardPosExpect)*ones(1,nSims+1), 'k-');
xlabel('time step');
ylabel('$$E\{|X_k|-|\hat{X}_k|\}$$','Interpreter','Latex');
title('Expected Absolute Error on Number of targets VS. Time step');
set(gca, 'Box', 'on');

figure(6);
subplot(211);whitebg('white'); cla; 
XCardAvg = sum(round(XCardHat) - XCardTrue, 2)'/nMCs;
XCardRms = sqrt(sum((round(XCardHat) - XCardTrue).^2, 2)'/nMCs-XCardAvg.^2);
plot(2:nSims, XCardAvg(2:end), 'k-','linewidth',2.5);
hold on;
plot(0:nSims, mean(XCardAvg)*ones(1,nSims+1), 'k-');
xlabel('time step');
ylabel('$$E\{\hat{n}_X-n_X\}$$','Interpreter','Latex');
title('Average Error of Target RFS Cardinality VS. Time step');
set(gca, 'Box', 'on');
subplot(212);whitebg('white'); cla
plot(2:nSims, XCardRms(2:end), 'k-','linewidth',2.5);
hold on;
plot(0:nSims, mean(XCardRms)*ones(1,nSims+1), 'k-');
xlabel('time step');
ylabel('RMSE');
title('Root Mean Square Error of Target RFS Cardinality VS. Time step');
set(gca, 'Box', 'on');

figure(7);
whitebg('white'); cla; hold on;
plot(2:nSims, XCardTrue(2:end,m), 'k-','linewidth',2.5);
XCardMean = mean(XCardHat,2)';
XCardVar = std(XCardHat, 0, 2)';
plot(2:nSims, XCardMean(2:end), 'ko','linewidth',2.5);
plot(2:nSims, XCardMean(2:end)-XCardVar(2:end), 'k--','linewidth',2.5);
plot(2:nSims, XCardMean(2:end)+XCardVar(2:end), 'k--','linewidth',2.5);
title('Cardinality mean VS time');
legend('true', 'mean', 'StDev');

figure(8);
whitebg('white'); cla; 
plot(2:nSims, mean(OSPA(2:end,:)'), 'k-d', 'linewidth', 2); 
hold on;
plot(2:nSims, mean(Hausdorf(2:end,:)'), 'k-+', 'linewidth', 1); 
plot(2:nSims, mean(OMAT(2:end,:,1)'), 'k-o', 'linewidth', 1); %#ok<*UDIM>
plot(2:nSims, mean(OMAT(2:end,:,2)'), 'k-*', 'linewidth', 1);
xlabel('time step');
ylabel('Metric');
title('Metric between true state RFS and Filter RFS VS. Time step');
legend('OSPA','Hausdorf','OMAT(p=1)', 'OMAT(p=2)');