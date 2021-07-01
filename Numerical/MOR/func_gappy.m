% ------------------------------------------------------------------------------
% Gappy function evaluation.
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function g = func_gappy(p,y,mor_object)
    out = p(y);
    
    f = dot(repmat(out,1,mor_object.g),mor_object.U_gap_S).';
    b = mor_object.M*f;
    g_tilde = sum(b'.*mor_object.U_gap,2);    
    g = zeros(size(y,1),1);
    g(mor_object.S) = out;
    g(mor_object.idx) = g_tilde(mor_object.idx);    
end
