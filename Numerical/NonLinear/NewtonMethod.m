% ------------------------------------------------------------------------------
% Newton method solver for nonlinear system.
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function [out, error] = NewtonMethod(F,x_l,J,tol)
    flag = 1; l = 1;
    while flag
        F_val = F(x_l);
        dx = J\(-F_val);
        x_l_1 = x_l + dx;
        error = norm(x_l_1 - x_l) + norm(F_val);
        if isnan(error)
            break;
        end    
        if error > tol && l < 30
            x_l = x_l_1;
            l = l+1;   
        else
            flag = 0;
        end
    end
    
    out = x_l;
end
