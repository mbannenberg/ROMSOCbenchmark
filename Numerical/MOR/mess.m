% ------------------------------------------------------------------------------
% Maximum Entropy Snapshot Sampling method.
%
% Copyright 2021 Fotios Kasolis (BUW, kasolis@uni-wuppertal.de)
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function [ Q ] = mess(X, tol)
if ( nargin == 1 )
  tol = 0.05;
end

n = size(X, 2);

% DISTANCE MATRIX
S = sum(X.^2, 1);
P = transpose(transpose(X)*X);
D = sqrt(abs(S + transpose(S) - 2.0*P));
D = D./max(D(:));

% RECURRENCE MATRIX
R = zeros(n,n);
R(D < tol) = 1;
R = tril(R);

% SELECTION PROCESS
Y = zeros(size(X));
itr = 1;
for cnt = 1:n
    idx  = find(R(:,cnt));
    if (isempty(idx))
        continue;
    else
        Y(:,itr) = mean(X(:,idx),2);
        R(:,idx) = 0.0;
        R(idx,:) = 0.0;
        itr = itr + 1;
    end
end
[Q,~] = qr(Y(:,1:itr-1),0);
