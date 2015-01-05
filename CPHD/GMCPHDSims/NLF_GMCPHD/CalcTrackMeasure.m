function [cpep, nabs] = CalcTrackMeasure(Track, XCLFilter, H, r)
% Routine used to calculate the track measure between trure track and filter
%  Tracks.
%
% Inputs:
%   Track -- Struct which save all the true track information.
%       Track.len --length of current track.
%       Track.t0 -- start time step index of current track.
%       Track.tf -- end time step index of current track.
%       Track.state -- state vector at current time step.
%       Track.id   -- string used to label current track.
%       Track.isactive -- flag indicates that current track is alive.
%       Track.offset -- offset temporary scalar in state stack vector.
%   XCLFilter -- Cell that stored all the extracted state from GMPHD.
%   H  -- Observation matrix.
%   r -- circular that used to calculate the CPEP.
%
% Outputs:
%   cpep -- Circular position error probability in every time step.
%   nabs -- absolute target number error.
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

if isempty(Track) || isempty(XCLFilter)
    return;
end
nMCs = size(Track(1).state,3);
nSims = size(Track(1).state, 2);
cpep = zeros(1, nSims);
nabs = zeros(1, nSims);
r2   = r^2;
for m=1:nMCs    
    for k=1:nSims
        X = [];
        for n=1:length(Track)
            if Track(n).len(m) && (k>=Track(n).t0(m)) && (k<=Track(n).tf(m))
                X = [X Track(n).state(:,k,m)]; %#ok<*AGROW>
            end
        end
    
        Xhat = XCLFilter{k,m};
        
        cp = 0;
        nX = size(X,2);
        nXhat = size(Xhat,2);
        if nX && nXhat % True RFS and Estimated RFS are both not empty set.
            for ii=1:nX
                Dis = sum((repmat(H*X(:,ii),1,nXhat)-H*Xhat).^2);
                if all(Dis>r2)
                    cp = cp + 1;
                end
            end
            cpep(k) = cpep(k) + cp/nX;
        elseif ~nX && ~nXhat  % True RFS and Estimated RFS are both empty set. 
                              % So this is correct estimation. CPEP stays.
            cpep(k) = cpep(k);
        else  % If any of X and Xhat is not empty, then miss estimation, 
              % cpep at current time is set to 1.
            cpep(k) = cpep(k) + 1;
        end
               
        nabs(k) = nabs(k) + abs(nX - nXhat);       
    
    end
end
nabs = nabs/nMCs;
cpep = cpep/nMCs;

