function [Combs, CurrentIndex] = GenCombinationIndex(p, n, m, RetainIndex, CurrentIndex, Combs)
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

if ~m
    Combs = [Combs, sort(CurrentIndex)];
    return;
end


% Choose the first element of Retain index as the first element of CurrentIndex.
if ~isempty(RetainIndex)
    CurrentIndex(p-m+1) = RetainIndex(1);
    RetainIndex = RetainIndex(2:end);
    [Combs, CurrentIndex] =...
        GenCombinationIndex(p, n-1, m-1, RetainIndex, CurrentIndex, Combs);

    CurrentIndex(1:p-m) = CurrentIndex(1:p-m);
    [Combs, CurrentIndex] = ...
        GenCombinationIndex(p, n-1, m,   RetainIndex, CurrentIndex, Combs);

end