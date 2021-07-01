% ------------------------------------------------------------------------------
% Reduced order multirate integration scheme using the BDF method.
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function [t, y] = ROMRBDF(E,A,func_p,func_r,y0,t,k,tol,m,mor_object)
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
    func_r_F = str2func(strcat('@(t) [', B{idx_F}, ']'));
    
    func_r = str2func(strcat('@(t) [', B{:}, ']'));
    
    B = strcat(regexprep(cellfun(@func2str, func_p, 'uni', 0), '^@\(x\)', ''), ';');
    func_p =  str2func(strcat('@(x) [', B{:}, ']'));
    func_p_F = str2func(strcat('@(x) [', B{idx_F}, ']'));
    func_p_gappy_selected = str2func(strcat('@(x) [', B{mor_object.S}, ']'));

    f = @(x_p,x,t) A*x + E*x_p + func_p(x) + func_r(t);
   
    f_F = @(x_p,x,t) A(idx_F,:)*x + E(idx_F,:)*x_p + func_p_F(x) + func_r_F(t);
        
    Joptions.diffvar = 1;
    Joptions.vectvars = [];
    Joptions.thresh = 1e-5*ones(M,1);
    Joptions.fac = [];

    f0 = func_p(y0);
    [J_p,Joptions.fac] = odenumjac(@(x) func_p(x), {y0}, f0, Joptions); 
    H = t(m+1) - t(1);
    J = A + 1/H*E + J_p;
    
    Joptions_F.diffvar = 1;
    Joptions_F.vectvars = [];
    Joptions_F.thresh = 1e-5*ones(size(idx_F,1),1);
    Joptions_F.fac = [];

    for n = m+1:m:N
        t_l = t(n);
        x_l = y(:,n-m);
        H = t(n) - t(n-m);
        f_BDF = @(x) f((x-x_l)/H,x,t_l);
        
        [x_l_1, error] = fsolveJac_full(f_BDF,x_l,J,tol,mor_object);
        
        if error > 1e-8 || isnan(error)            
            f0 = func_gappy(func_p_gappy_selected,x_l,mor_object);
            [J_p,Joptions.fac] = odenumjac(@(x) func_gappy(func_p_gappy_selected,x,mor_object), {x_l}, f0, Joptions); 
            J = A+1/H*E+J_p;
            [x_l_1, error] = NewtonMethod(f_BDF,x_l,J,tol); %NR(F,x_l,x0,h,t_l,J,tol);
            if isnan(error)
                a = 3;
            end
        end

        if m == 1
            y(:,n) = x_l_1;
        else
            y(idx_S,n+1-m:n) = repmat(x_l_1(idx_S),1,m);

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


function [out, error, k] = fsolveJac_full(f,x,J,tol,mor_object)
    flag = 0;
    k = 0;
    
    w0 = x; w_r_k = zeros(size(mor_object.U_r'*x));
    while flag == 0
        temp = w0+mor_object.U_r*w_r_k;
        F_k = f(temp);
        w_r_k_1 = w_r_k + (J*mor_object.U_r)\(-F_k);
        error = norm(w_r_k_1 - w_r_k);
        
        if error > tol && k < 100
            w_r_k = w_r_k_1;
            k = k+1; 
        else
            flag = 1;
            if error > tol || isnan(error)
            end
        end
    end 
    out = w0 + mor_object.U_r*w_r_k_1;
end