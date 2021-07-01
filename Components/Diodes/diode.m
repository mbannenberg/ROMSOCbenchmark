% -----------------------------------------------------------------------------
% Shockley diode characteristic funtion
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function out = diode(I_s,v1,v2)
      out = I_s*(exp((v1-v2)/0.0256)-1);
end