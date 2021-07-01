% ------------------------------------------------------------------------------
% Gappy proper orthogonal decomposition method.
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function [ S, M_ij, U_gap_S, Idx]  = gappypod(U)
    [n,m] = size(U) ;
    [~,~,P] = qr(U','vector'); 
    S = P(1:m); 
    U_gap_S = U(S,:);

    n_k = zeros(n,1);
    n_k(S) = 1;
    Idx = find(n_k == 0);

    M_ij = zeros(m,m);
    for i = 1:m
        for j = 1:m
            M_ij(i,j) = norm_n(n_k,U(:,i),U(:,j));
        end
    end
    M_ij = inv(M_ij);
end

function out = norm_n(n,u,v)
    out = dot(n.*u,n.*v);
end