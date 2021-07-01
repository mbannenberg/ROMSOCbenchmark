% ------------------------------------------------------------------------------
% Voltage source
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function out = v_in_delay(t)
    if (t < 15e-9)
        out = 0;
        return;
    end

    t = t*1e7;
    multi = 10*1/2*180;
    phase = 0.1;
    t = mod(t,phase*1.5);
    u_d_p = 0.02;
    if(t < u_d_p*phase)
        out = multi*t;
    elseif (t <= (1- u_d_p)*phase)
        out = multi*u_d_p*phase;
    elseif (t <= phase)
        out = (phase - t)*multi;
    else
        out = 0;
    end
end