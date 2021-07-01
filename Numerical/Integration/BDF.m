% ------------------------------------------------------------------------------
% BDF method integration scheme.
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function [t, y] = BDF(E,A,func_p,func_r,y0,t,k,tol)
    N = size(t,2);
    h = t(2) - t(1);
    y = zeros(size(y0,1),size(t,2));
    y(:,1) = y0;

    B = strcat(regexprep(cellfun(@func2str, func_p, 'uni', 0), '^@\(x\)', ''), ';');
    func_p =  str2func(strcat('@(x) [', B{:}, ']'));
    
    B = strcat(regexprep(cellfun(@func2str, func_r, 'uni', 0), '^@\(t\)', ''), ';');
    func_r =  str2func(strcat('@(t) [', B{:}, ']'));
    
    f = @(x_p,x,t) A*x + E*x_p + func_p(x) + func_r(t);
        
    Joptions.diffvar = 1;
    Joptions.vectvars = [];
    Joptions.thresh = 1e-5*ones(size(y0));
    Joptions.fac = [];

    f0 = func_p(y0);
    [J_p,Joptions.fac,nF] = odenumjac(@(x) func_p(x), {y0}, f0, Joptions); 
    J = A+1/h*E+J_p;
    J_base = A+1/h*E;


    for i = 2: N
        t_l = t(i);
        x_l = y(:,i-1);
        
        % Define derivative approximation
        f_BDF = @(x) f((x-x_l)/h,x,t_l);
        [x_l_1, error] = NewtonMethod(f_BDF,x_l,J,tol);


        if error > 1e-8 || isnan(error)
            f0 = func_p(x_l);
            [J_p,Joptions.fac,nF] = odenumjac(@(x) func_p(x), {x_l}, f0, Joptions); 
            J = J_base+J_p;
            [x_l_1, error] = NewtonMethod(f_BDF,x_l,J,tol); 
        end

        if (isnan(rcond(J)) || rcond(J) < 1e-15)
            y = y(:,1:i-1);
            t = t(:,1:i-1);
            disp('Error in Jac');
            break;
        end
        y(:,i) = x_l_1;
    end
end