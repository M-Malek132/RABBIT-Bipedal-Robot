%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE VARIABLE: relative stance-foot -> hip vector angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function th = theta_of_q(q)
    hip  = q(1:2);
    foot = P_st(q);
    rel  = hip - foot;
    th   = atan2(rel(1), rel(2));
end