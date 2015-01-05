function dist= omat_dist(X,Y,p)
%
%Compute Optimal Mass Transfer Metric between two finite sets X and Y
%as described in the reference
%[1] D. Schuhmacher, B.-T. Vo, and B.-N. Vo, "A consistent metric for performance evaluation in multi-object 
%filtering," IEEE Trans. Signal Processing, Vol. 56, No. 8 Part 1, pp. 3447?3457, 2008.
%
%Inputs: X,Y--   matrices of column vectors
%        p  --   p-parameter for the metric.
%Output: scalar OMAT distance between X and Y        
%Note: the Euclidean 2-norm is used as the "base" distance on the region
%

if isempty(X) && isempty(Y)
    dist = 0;
    return;
end

if isempty(X) || isempty(Y)
    if isempty(X)  % If any set between X and Y is empty, 
                   % return the mean distance of the noempty set.
        dist = sum(sqrt(sum(Y.^2)));
    end
    if isempty(Y)
        dist = sum(sqrt(sum(X.^2)));
    end
   
    return;
end


%Calculate sizes of the input point patterns
n = size(X,2);
m = size(Y,2);
d = gcd(n,m);
n1 = m/d;
m1 = n/d;
if m~=n
    XX = [];
    YY = [];
    for ii=1:n
        XX = [XX repmat(X(:,ii), [1 n1])];
    end
    XX = repmat(XX, [1 n*n1]);
    
    for jj=1:m
        YY = [YY repmat(Y(:,jj), [1, m1])]; %#ok<*AGROW>
    end
    YY = reshape(repmat(YY,[n*n1 1]),[size(Y,1) n*m*n1*m1]);
else
    %Calculate cost/weight matrix for pairings - fast method with vectorization
    XX = repmat(X,[1 m]);
    YY = reshape(repmat(Y,[n 1]),[size(Y,1) n*m]);
end


if ~isinf(p)
    D  = reshape(sqrt(sum((XX-YY).^2)).^p,[n*n1 m*m1]);    
else
    D  = reshape(sqrt(sum((XX-YY).^2)),[n*n1 m*m1]);    
end

% %Calculate cost/weight matrix for pairings - slow method with for loop
% D= zeros(n,m);
% for j=1:m
%     D(:,j)= sqrt(sum( ( repmat(Y(:,j),[1 n])- X ).^2 )');
% end
% D= min(c,D).^p;

%Compute optimal assignment and cost using the Hungarian algorithm
[assignment,cost]= Hungarian(D); %#ok<*ASGLU>

%Calculate final distance
dist = (cost/n/n1).^(1/p);



    