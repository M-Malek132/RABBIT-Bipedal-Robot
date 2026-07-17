function g = dtheta_dq_of(q)
    hip  = q(1:2);
    foot = P_st(q);
    rel  = hip - foot;
    r2   = max(rel(1)^2 + rel(2)^2, 1e-10);   % avoid divide-by-zero

    dhip_dq = [1 0 zeros(1,5); 0 1 zeros(1,5)];
    Jst     = J_st(q);
    drel_dq = dhip_dq - Jst;

    g = ( rel(2)*drel_dq(1,:) - rel(1)*drel_dq(2,:) ) / r2;
end