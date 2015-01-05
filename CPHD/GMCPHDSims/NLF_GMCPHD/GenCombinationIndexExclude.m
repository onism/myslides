function [Combs, CurrentIndex, IndexNext] = GenCombinationIndexExclude(Combs, CurrentIndex, IndexRetain, lvl)
% Function used to calculate all the Combination of Indexing series.
%
% Inputs:
%     Combs -- Current Combinations collections.
%     CurrentIndex  -- Current Combination index.
%     IndexRetain -- Index restain.(vector)
%     lvl -- Level the iterative routine has to run.(scalar, nonegative integer)
% Notice: The value of lvl must be less than or equal to  length(IndexRestain).
%     For speed computation purpose, we will NOT validate that. 
%     You should set the initial value in illegal value.
%
% Outputs:
%    Combs -- Current Combinations collections after current run.
%    CurrentIndex -- Current Combination Index.
%    IndexNext -- Index for next run.
%
%
%   Copyright 2012-2013 N.U.D.T. 
%   Author:   David Lieu. 
%   $Revision: 1.0 $  
%   $Date: 2012/04/03 08:53:21 $
%
N = length(CurrentIndex);
if ~lvl
    if isempty(Combs)
        Combs = sort(CurrentIndex);
        return;
    end
    
    flag = all(sum((Combs-repmat(sort(CurrentIndex),[1,size(Combs,2)])).^2));
    if flag
        Combs = [Combs, sort(CurrentIndex)];
    end
    return;    
end

for n=1:length(IndexRetain)
    CurrentIndex(N-lvl+1) = IndexRetain(n); %#ok<*AGROW>
    IndexNext = [IndexRetain(1:n-1);IndexRetain(n+1:end)];
    [Combs, CurrentIndex] = ...
        GenCombinationIndexExclude(Combs, CurrentIndex, IndexNext, lvl-1);
end

