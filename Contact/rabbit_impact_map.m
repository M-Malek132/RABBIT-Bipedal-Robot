function x_plus = rabbit_impact_map(x_minus)
    q  = x_minus(1:7);
    dq = x_minus(8:14);

    D = M(q);
    J = J_sw(q); 

    % System of equations: 
    % [ M  -J' ] [ dq_plus ] = [ M * dq_minus ]
    % [ J   0  ] [ lambda  ]   [ 0            ]
    
    A = [D, -J';
         J,  zeros(size(J,1))];
    b = [D*dq;
         zeros(size(J,1),1)];

    sol = A\b;
    dq_plus = sol(1:7);

    % Note: q does not change during impact, only dq
    x_plus = [q; dq_plus];
end
