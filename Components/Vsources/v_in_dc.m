% ------------------------------------------------------------------------------
% Voltage source
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function out = v_in_dc(t,val)
    t = t*1e9;
    if (t < 1)
        out = val * t;
    else
        out = val;
    end
end