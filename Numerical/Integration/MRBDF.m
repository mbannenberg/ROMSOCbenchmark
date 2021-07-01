% ------------------------------------------------------------------------------
% Multirate integration scheme using BDF method.
% 
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function [t, y] = MRBDF(E,A,func_p,func_r,y0,t,k,tol,m)
    N = size(t,2); M = size(y0,1);
    
    y = zeros(M,N); y(:,1) = y0;
    
    idx_F = [1:3 M]';
    idx_S = [4:M-1]';
        
    assert(size(idx_F,1)+size(idx_S,1) == M);
    
    B_F = zeros(size(idx_F,1),M); B_S = zeros(size(idx_S,1),M);
    
    for i = 1:size(idx_F,1)
        B_F(i,idx_F(i)) = 1;
    end
    for i = 1:size(idx_S,1)
        B_S(i,idx_S(i)) = 1;
    end
    
    % Make partitioned functions for each partition.
    B = strcat(regexprep(cellfun(@func2str, func_r, 'uni', 0), '^@\(t\)', ''), ';');
    
    func_r = str2func(strcat('@(t) [', B{:}, ']'));
    
    B = strcat(regexprep(cellfun(@func2str, func_p, 'uni', 0), '^@\(x\)', ''), ';');
    func_p =  str2func(strcat('@(x) [', B{:}, ']'));    
    
    % Construct the implicit MNA function of the circuit.
    f = @(x_p,x,t) A*x + E*x_p + func_p(x) + func_r(t);

    % Define jacobian calculation parameters and compute first Jacobian.
    Joptions.diffvar = 1;
    Joptions.vectvars = [];
    Joptions.thresh = 1e-5*ones(M,1);
    Joptions.fac = [];
    f0 = func_p(y0);
    [J_p,Joptions.fac] = odenumjac(@(x) func_p(x), {y0}, f0, Joptions); 
    H = t(m+1) - t(1);
    J = A+1/H*E+J_p;
    
    Joptions_F.diffvar = 1;
    Joptions_F.vectvars = [];
    Joptions_F.thresh = 1e-5*ones(size(idx_F,1),1);
    Joptions_F.fac = [];

    for n = m+1:m:N
        t_l = t(n);
        x_l = y(:,n-m);
        H = t(n) - t(n-m);
        
        f_BDF = @(x) f((x-x_l)/H,x,t_l);
        
        [x_l_1, error] = NewtonMethod(f_BDF,x_l,J,tol);
         
        if error > 1e-8 || isnan(error)
            f0 = func_p(x_l);
            [J_p,Joptions.fac,nF] = odenumjac(@(x) func_p(x), {x_l}, f0, Joptions); 
            J = A+1/H*E+J_p;
            [x_l_1, error] = NewtonMethod(f_BDF,x_l,J,tol);
        end
        y(idx_S,n) = x_l_1(idx_S);

        if m == 1
            y(idx_F,n) = x_l_1(idx_F);
        else
            y(idx_S,n+1-m:n-1) = repmat(x_l_1(idx_S),1,m-1);
            for l = 0:m-1
                
                h = t(n-m+l+1) - t(n - m + l);
                t_l_1 = t(n-m+l+1);
                
                x_F_l = y(idx_F,n-m+l);
                x_S_l = y(idx_S,n-m+l);
                x_S_l_1 = y(idx_S,n-m+l+1);

                f_F_BDF = @(x) f_F((B_F'*(x-x_F_l)+ B_S'*(x_S_l_1-x_S_l))/h, B_F'*x + B_S'*x_S_l_1, t_l_1);
                
                f0_F = f_F_BDF(x_F_l);
                                
                [J_F,Joptions_F.fac] = odenumjac(@(x) f_F_BDF(x), {x_F_l}, f0_F, Joptions_F); 

                [x_F_l_1, error] = NewtonMethod_F(f_F_BDF,x_F_l,J_F,tol);
                                
                y(idx_F,n-m+l+1) = x_F_l_1;                
            end
         end
    end
end