function esv = ElementarySymmetricFunc(Z, p)
%ELEMENTARYSYMMETRICFUNC is the elementary symmetric function.
%
%  Inputs:
%    Z --- the list of data. It's a column vector. If it's a matrix,
%          each column will be treated as a list of data.
%    p -- Order of elementary symmetric function.
%
%  Ouputs:
%    esv-- Elementary symmetric value.
%
%   Copyright 2012-2013 N.U.D.T. 
%   Author:   David Lieu. 
%   $Revision: 1.0 $  
%   $Date: 2012/03/20 16:42:24 $
%
%
[M, N] = size(Z);
esv = zeros(1, N);
CombsIndex = GenCombinationIndex(p,M,p,(1:M)',zeros(p,1),[]);
for n=1:size(CombsIndex,2)
    Idx = CombsIndex(:,n);
    esv = esv + prod(Z(Idx,:),1);
end





