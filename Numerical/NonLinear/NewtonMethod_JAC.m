% ------------------------------------------------------------------------------
% Newton method solver for nonlinear system. Jacobian recomputed each step.
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function [out, error] = NewtonMethod_JAC(F,x_l,J,tol)
    flag = 1; l = 1;
    Joptions.diffvar = 1;
    Joptions.vectvars = [];
    Joptions.thresh = 1e-13*ones(size(x_l,1),1);
    Joptions.fac = [];
    while flag
        F_val = F(x_l);
        [J,Joptions.fac] = odenumjac(F, {x_l}, F_val, Joptions); 

        dx = J\(-F_val);
        x_l_1 = x_l + dx;
        error = norm(x_l_1 - x_l) + norm(F_val);
        if isnan(error)
            break;
        end    
        if error > tol && l < 300
            x_l = x_l_1;
            l = l+1;   
        else
            flag = 0;
        end
    end
    
    out = x_l;
end
