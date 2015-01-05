function esv = ElemSymmFunc(Z, p)
%ElemSymmFunc is the elementary symmetric function.
%
%  Inputs:
%    Z -- the list of data. It's a column vector. If it's a matrix,
%          each column will be treated as a list of data.
%    p -- Order of elementary symmetric function.
%
%  Ouputs:
%    esv-- Elementary symmetric value.
%
%   Copyright 2012-2013 N.U.D.T. 
%   Author:   David Lieu. 
%   $Revision: 1.0 $.  
%   $Date: 2012/03/20 16:42:24 $
%
%
[M,N] = size(Z);
X     = zeros(M+1,N, class(Z));
for n=1:N
    X(:,n) = poly(Z(:,n));
end
X = X./repmat(X(1,:),[M+1,1]);

I = (-1).^(0:M)';
ESV = X.*I;

if nargin==1    
    esv = ESV';    
else
    esv   = zeros(size(p));
    %if all(p>=0 & p<=M)
    esv(p>=0&p<=M) = ESV(p(p>=0&p<=M)+1,:);
    %end
end 
