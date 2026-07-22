function [x_plus, impulse] = rabbit_impact_map(x_minus)
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

    % lambda here is the CONTACT IMPULSE [Ns] at the swing foot (world [x;z]).
    % Returned as an optional 2nd output; existing callers use only x_plus.
    impulse = sol(8:end);

    % Note: q does not change during impact, only dq
    x_plus = [q; dq_plus];
end
