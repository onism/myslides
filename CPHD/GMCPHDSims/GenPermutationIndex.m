function [Perms, CurrentIndex, IndexNext] = GenPermutationIndex(Perms, CurrentIndex, IndexRetain, lvl)
% Function used to calculate all the permutation of Indexing series.
%
% Inputs:
%     Perms -- Current permutations collections.
%     CurrentIndex  -- Current Permutation index.
%     IndexRetain -- Index restain.(vector)
%     lvl -- Level the iterative routine has to run.(scalar, nonegative integer)
% Notice: The value of lvl must be less than or equal to  length(IndexRestain).
%     For speed computation purpose, we will NOT validate that. 
%     You should set the initial value in illegal value.
%
% Outputs:
%    Perms -- Current permutations collections after current run.
%    CurrentIndex -- Current Permutation Index.
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
    Perms = [Perms, CurrentIndex];
    return;
end

for n=1:length(IndexRetain)
    CurrentIndex(N-lvl+1) = IndexRetain(n); %#ok<*AGROW>
    IndexNext = [IndexRetain(1:n-1);IndexRetain(n+1:end)];
    [Perms, CurrentIndex] = ...
        GenPermutationIndex(Perms, CurrentIndex, IndexNext, lvl-1);
end

