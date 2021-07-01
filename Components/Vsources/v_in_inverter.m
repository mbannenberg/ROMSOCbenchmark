% ------------------------------------------------------------------------------
% Voltage source
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function out = v_in_inverter(t)
    t = t*1e7;
    multi = 7200;
    t = mod(t,0.002);
    if(t < 0.00025)
        out = multi*t;
    elseif (t <= 0.00075)
        out = multi*0.00025;
    elseif (t <= 0.001)
        out = (0.001 - t)*multi*0.00025/0.00025;
    else
        out = 0;
    end
end