%% Initialization of True Track with STRUCT structure.
% Setting state segment of initial Track.
SX(:,1)  = [-800 20 -200 -5]';  
SX(:,2)  = [-800 12.5 -200 7]';
SX(:,3)  = [-800 2.5 -200 10.5]';
SX(:,4)  = [-200 16 800 -9.7]';
SX(:,5)  = [-200 -2.5 800 -14.6]';
SX(:,6)  = [-200 17.5 800 -5]';
SX(:,7)  = [0    0 0 -10]';
SX(:,8)  = [0    7.5 0 -5]';
SX(:,9)  = [0    -20 0 -15]';
SX(:,10) = [400  -10.5 -600 5]';
SX(:,11) = [400  -2.5  -600 10.2]';
SX(:,12) = [400  -7.5  -600 -4.5]';

% Setting Life span segment of Initial Track.
SP(:,1)  = [1   70]';    
SP(:,2)  = [40 100]';
SP(:,3)  = [60 100]';
SP(:,4)  = [40 100]';
SP(:,5)  = [60 100]';
SP(:,6)  = [80 100]';
SP(:,7)  = [1   70]';
SP(:,8)  = [20 100]';
SP(:,9)  = [80 100]';
SP(:,10) = [1  100]';
SP(:,11) = [20 100]';
SP(:,12) = [20 100]';

% Setting type segment of initial track.
ST(:,1)  = TYPE_BIRTH;
ST(:,2)  = TYPE_BIRTH;
ST(:,3)  = TYPE_BIRTH;
ST(:,4)  = TYPE_BIRTH;
ST(:,5)  = TYPE_BIRTH;
ST(:,6)  = TYPE_BIRTH;
ST(:,7)  = TYPE_BIRTH;
ST(:,8)  = TYPE_BIRTH;
ST(:,9)  = TYPE_BIRTH;
ST(:,10) = TYPE_BIRTH;
ST(:,11) = TYPE_BIRTH;
ST(:,12) = TYPE_BIRTH;

% Setting isactive segment of initial track.
SA     = zeros(1,12);
SA(1)  = 1;
SA(7)  = 1;
SA(10) = 1;

% Setting the offset segment of initial track.
SO(1)  = 1;
SO(2)  = 1;
SO(3)  = 1;
SO(4)  = 2;
SO(5)  = 2;
SO(6)  = 2;
SO(7)  = 3;
SO(8)  = 3;
SO(9)  = 3;
SO(10) = 4;
SO(11) = 4;
SO(12) = 4;

%% Track Initialization
for nt=1:12
    Track(nt).len      = zeros(1, nMCs);  %#ok<*SAGROW> % length of current track.    
    Track(nt).t0       = SP(1,nt)*ones(1, nMCs);    % start time step index of track.
    Track(nt).tf       = SP(2,nt)*ones(1, nMCs);  % final time step index of track.
    Track(nt).state    = repmat([zeros(nDimX, SP(1,nt)-1) SX(:,nt) zeros(nDimX, nSims-SP(1,nt))],[1,1,nMCs]);  % State vector of current track in different Monte Carlo runs.    
    Track(nt).model    = MODEL_CV*ones(1, nSims);   % model type that current track follows.  (MODEL_CV, MODEL_CT, MODEL_CA)
    if nt<10
        Track(nt).id   = ['BT0' num2str(nt)];       % Track identity or label.
    else
        Track(nt).id   = ['BT' num2str(nt)];
    end
    Track(nt).type     = ST(:,nt); % Type can be set to TYPE_BIRTH or TYPE_SPAWN.
    Track(nt).isactive = SA(nt)*ones(1, nMCs);      % isactive is the flag whether current track is being updated.
    Track(nt).offset   = SO(nt);
end
clear SX SP ST SA SO;
