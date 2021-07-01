% ------------------------------------------------------------------------------
% Epsilon estimation procedure.
%
% Copyright 2021 Fotios Kasolis (BUW, kasolis@uni-wuppertal.de 
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function out = epsilon_procedure(X)
    n = size(X,2);

    S = sum(X.^2, 1);
    P = transpose(transpose(X)*X);
    D = sqrt(abs(S + transpose(S) - 2.0*P));
    D = D./max(D(:));
    
    epsilon = linspace(0.005,0.99,300);

    for cnt = 1:numel(epsilon)
        tol = epsilon(cnt);
        R = zeros(n,n);
        R(D < tol) = 1;
        C(cnt) = sum(R(:))./(n*n);
    end

    ln_P_E_over_ln_eps = log(C)./log(epsilon);

    x = epsilon;
    r=randdf([30000],[C;x],'cdf'); % generate random numbers
    r = r(r>=0);
    avg_ln_y = mean(-log(r)).^(-1);
    [a, idx] = min(abs(avg_ln_y - ln_P_E_over_ln_eps));
    out = epsilon(idx);
end