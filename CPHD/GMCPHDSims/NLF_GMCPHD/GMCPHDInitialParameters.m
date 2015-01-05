%% Scripts used to initialize the current tracking view.
%  MATLAB_DEBUG = 1;  % When you just want to run the GMPHD_LinearGaussianMotionSim just once.
                      % You can just set this variable to 1. When you want to varify some
                      % parameter to obtain the perfornamce factors you can set it to 0 explicitly.
strMarker = 'px+*sdhv^><o';
flagSpecified = 1;
flagShowWaitbar = 0;

%% Number of simulation steps and Monte Carlo runs.
T = 1;        %Sampling interval.
nSims = 100;  % Number of simulation step.
nMCs  = 64;    % number of Monte Carlo runs.

%% Seting basic dimension of state and observation vector. 
nDimX = 4;    % Dimension of State vector.
nDimZ = 2;    % Dimension of Measurement vector.

%% Surveillance region.
xmin = -1000;  % Minimum position in x axis.
xmax = 1000;   % Maximum position in x axis.
ymin = -1000;  % Minimum position in y axis.
ymax = 1000;   % Maximum position in y axis.
vol = (xmax-xmin)*(ymax-ymin);  % Volume of surveillance region.
u = 1/vol;     % uniform distribution of clutter.

%% Clutter density.
lambdac = 12.5e-6;%1e-6;%10e-7;%Poisson clutter spatial distribution density.

%% CPEP parameter.
cpep_radius = 20;  % radius used to calculate the CPEP.

%% process noise and observation noise setting.
sigmav = 15;        % standard variance of process noise.
sigmaw = 10;       % standard variance of measurement noise.
sigmau = deg2rad(1);  % standard variance of angle rate noise.

%% Survival probability and detection probability.
probsv = 1;      % Probability of target survival.
probdt = 0.98;   % Probability of detection.
pspar  = 0.99;   % Probability of target survival as a GM-PHD filter parameter.
pdpar  = 0.98;   % Probability of detection as a GM-PHD filter parameter.

%% Threshold and maximum number of gaussian components sets
thtrunc = 1e-5;  % threshold to truncate Gaussian component with small weight.
thunion = 4;     % threshold to union those Gaussian component with near means.
nGMMax = 100;    % Maximum Number of Gaussian component.
%kspawn = 66;     % Time step at which Spawning target appears.
nCardMax = 100;

%% System model and Observation model.
F = kron(eye(2),[1 T; 0 1]);
Q = sigmav*sigmav*kron(eye(2),[T^4/4 T^3/2; T^3/2 T^2]);
H = kron(eye(2),[1 0]);
R = sigmaw*sigmaw*eye(2);

alpha_man = 0;
w = pi/18;   % angle rate.

Fcv = F;
Qcv = Q;
Fct = inline(['[1, T*sinc(w*T/pi), 0, -T^2*w/2*(sinc(w*T/2/pi)^2), 0;' ...
              '0, cos(w*T), 0, -sin(w*T), 0;' ...
              '0, T^2*w/2*(sinc(w*T/2/pi)^2), 1, T*sinc(w*T/pi), 0;' ...
              '0, sin(w*T), 0, cos(w*T), 0;'...
              '0, 0, 0, 0, exp(-alpha*T)];'],...
              'w', 'T', 'alpha');  % Transition Matrix in Constant Turn Model.
Fjb = inline(['[1,T*sinc(w*T/pi),0,-T^2*w/2*(sinc(w*T/2/pi)^2),'...
              '1/w/w+(w*T-1)*cos(w*T)/w/w-(w*T+1)*sin(w*T)/w/w;' ...
              '0,cos(w*T),0,-sin(w*T),-T*cos(w*T)-T*sin(w*T);' ...
              '0,T^2*w/2*(sinc(w*T/2/pi)^2),1,T*sinc(w*T/pi),'...
              '-1/w/w+(w*T+1)*cos(w*T)/w/w+(w*T-1)*sin(w*T)/w/w;' ...
              '0,sin(w*T),0,cos(w*T),T*cos(w*T)-T*sin(w*T);'...
              '0,0,0,0,exp(-alpha*T)];'],...
              'w', 'T', 'alpha');  % Transition Matrix in Constant Turn Model.          
Qct = sigmaw*kron(eye(2), [T^2/2;T]);
Qct = Qct*Qct';   % Variance matrix of process noise in Constant Turn Model.
Qct = [Qct zeros(size(Qct,1),1);zeros(1,size(Qct,1)) sigmau^2];
%% Gaussian Mixture Model for birth target.
GMM_birth(1).omega = 0.03;
GMM_birth(1).mean = [-800 0 -200 0]';
GMM_birth(1).variance = diag([10 10 10 10]);

GMM_birth(2).omega = 0.03;
GMM_birth(2).mean = [-200 0 800 0]';
GMM_birth(2).variance = diag([10 10 10 10]);

GMM_birth(3).omega = 0.03;
GMM_birth(3).mean = [0 0 0 0]';
GMM_birth(3).variance = diag([10 10 10 10]);

GMM_birth(4).omega = 0.03;
GMM_birth(4).mean = [400 0 -600 0]';
GMM_birth(4).variance = diag([10 10 10 10]);
%% Gaussian Mixture Model for birth target.
% GMM_spawn(1).omega = 0.03;
% GMM_spawn(1).mean = [0 0 0 0]';
% GMM_spawn(1).variance = diag([100 400 100 400]);
% GMM_spawn(1).trans = F;
GMM_spawn = [];  % There is no spawning target in CPHD filter.

SetTargetInitialTracks;
CheckInitialConfig;
