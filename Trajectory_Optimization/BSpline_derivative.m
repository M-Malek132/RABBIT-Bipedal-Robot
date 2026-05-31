function [dN] = BSpline_derivative(n, p, u)
% Compute derivative of B-spline basis functions analytically
%
% d/ds N_{i,p}(s) = p * [ N_{i,p-1}(s) / (u_{i+p} - u_i) 
%                        - N_{i+1,p-1}(s) / (u_{i+p+1} - u_{i+1}) ]
%
% For boundaries (u=0 or u=1), uses numerical derivative to avoid
% 0/0 singularities from repeated knots in clamped knot vector.
%
% Inputs:
%   n: number of data points (control points = n+1)
%   p: degree of B-spline  
%   u: parameter value [0,1]
%
% Output:
%   dN: 1 x (n+1) vector of basis function derivatives

    if p == 0
        dN = zeros(1, n+1);
        return;
    end

    m = n + p + 1;
    
    U = zeros(1, m+1);
    U(:, p+1:n+2) = linspace(0, 1, m+1-2*p);
    U(:, n+2:m+1) = 1;
    
    % Handle boundaries numerically (repeated knots cause 0/0)
    if u <= U(1)
        ds = 1e-8;
        N0 = BSpline(n, p, 0);
        N1 = BSpline(n, p, ds);
        dN = (N1 - N0) / ds;
        return;
    end
    
    if u >= U(end)
        ds = 1e-8;
        N0 = BSpline(n, p, 1-ds);
        N1 = BSpline(n, p, 1);
        dN = (N1 - N0) / ds;
        return;
    end
    
    % Interior: 0 < u < 1 - use recurrence
    N_lower = zeros(1, n+2);
    
    % Degree 0
    for i = 0:n+1
        idx = i + 1;
        if idx+1 <= length(U)
            if u >= U(idx) && u < U(idx+1)
                N_lower(idx) = 1;
                break;
            end
        end
    end
    
    % Recurse to degree p-1
    for d = 1:p-1
        N_temp = zeros(1, n+2);
        for i = 0:(n+2-d)
            idx = i + 1;
            if idx > length(N_lower), break; end
            
            term1 = 0; term2 = 0;
            
            if N_lower(idx) ~= 0 && i+d+1 <= length(U)
                denom1 = U(i+d+1) - U(i+1);
                if denom1 ~= 0
                    term1 = ((u - U(i+1)) / denom1) * N_lower(idx);
                end
            end
            
            if idx+1 <= length(N_lower) && N_lower(idx+1) ~= 0
                if i+d+2 <= length(U) && i+2 <= length(U)
                    denom2 = U(i+d+2) - U(i+2);
                    if denom2 ~= 0
                        term2 = ((U(i+d+2) - u) / denom2) * N_lower(idx+1);
                    end
                end
            end
            
            N_temp(idx) = term1 + term2;
        end
        N_lower = N_temp;
    end
    
    % Derivative formula
    dN = zeros(1, n+1);
    for i = 0:n
        idx = i + 1;
        term1 = 0; term2 = 0;
        
        if N_lower(idx) ~= 0 && i+p+1 <= length(U)
            denom1 = U(i+p+1) - U(i+1);
            if denom1 ~= 0
                term1 = p * N_lower(idx) / denom1;
            end
        end
        
        if idx+1 <= length(N_lower) && N_lower(idx+1) ~= 0
            if i+p+2 <= length(U) && i+2 <= length(U)
                denom2 = U(i+p+2) - U(i+2);
                if denom2 ~= 0
                    term2 = p * N_lower(idx+1) / denom2;
                end
            end
        end
        
        dN(idx) = term1 - term2;
    end
    
end