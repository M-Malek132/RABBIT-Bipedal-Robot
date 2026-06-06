% test_bspline_derivative.m
% Debug and verify BSpline_derivative against numerical differentiation

clear; clc; close all;

%% Test Case 1: Simple case with known derivative
fprintf('=== TEST CASE 1: Simple cubic B-spline ===\n');

n = 5;  % 6 control points
p = 3;  % cubic

% Test at several points
s_test = [0, 0.25, 0.5, 0.75, 1.0];
ds = 1e-7;

for s = s_test
    fprintf('\nAt s = %.4f:\n', s);
    
    % Get analytical derivative
    try
        dN = BSpline_derivative(n, p, s);
        fprintf('  Analytical dN: [');
        fprintf('%.4f ', dN);
        fprintf(']\n');
    catch ME
        fprintf('  ERROR in analytical: %s\n', ME.message);
    end
    
    % Numerical derivative for comparison
    if s <= 0
        s0 = s;
        s1 = s + ds;
    elseif s >= 1
        s0 = s - ds;
        s1 = s;
    else
        s0 = s - ds/2;
        s1 = s + ds/2;
    end
    
    N0 = BSpline(n, p, s0);
    N1 = BSpline(n, p, s1);
    dN_num = (N1 - N0) / (s1 - s0);
    
    fprintf('  Numerical dN:  [');
    fprintf('%.4f ', dN_num);
    fprintf(']\n');
    
    % Check sum property: sum of derivatives should be 0
    if exist('dN', 'var')
        fprintf('  Sum analytical: %.10f (should be ~0)\n', sum(dN));
        fprintf('  Sum numerical:  %.10f (should be ~0)\n', sum(dN_num));
        
        error_val = max(abs(dN - dN_num));
        fprintf('  Max error: %.2e\n', error_val);
    end
end

%% Test Case 2: Verify against known B-spline properties
fprintf('\n=== TEST CASE 2: B-spline derivative properties ===\n');

n = 7;  % 8 control points  
p = 3;  % cubic

% Property: sum of derivative basis functions should be 0
s_check = linspace(0, 1, 100);
sum_derivatives = zeros(size(s_check));

for i = 1:length(s_check)
    dN = BSpline_derivative(n, p, s_check(i));
    sum_derivatives(i) = sum(dN);
end

fprintf('Maximum absolute sum of derivatives: %.2e (should be near 0)\n', ...
        max(abs(sum_derivatives)));

if max(abs(sum_derivatives)) > 1e-10
    fprintf('WARNING: Sum of derivatives not zero!\n');
end

%% Test Case 3: Visual comparison with numerical derivative
fprintf('\n=== TEST CASE 3: Visual comparison ===\n');

n = 7;
p = 3;

% Create a specific set of control points
CP = [0.2, -0.5, 0.3, -0.1, 0.6, -0.3, 0.4, -0.2]';

% Fine grid for evaluation
s_vec = linspace(0, 1, 500);
analytical_deriv = zeros(size(s_vec));
numerical_deriv = zeros(size(s_vec));
ds_num = 1e-6;

for i = 1:length(s_vec)
    s = s_vec(i);
    
    % Analytical
    dN = BSpline_derivative(n, p, s);
    analytical_deriv(i) = dN * CP;
    
    % Numerical
    if s < ds_num
        s0 = s;
        s1 = s + ds_num;
    elseif s > 1 - ds_num
        s0 = s - ds_num;
        s1 = s;
    else
        s0 = s - ds_num/2;
        s1 = s + ds_num/2;
    end
    
    N0 = BSpline(n, p, s0);
    N1 = BSpline(n, p, s1);
    numerical_deriv(i) = (N1 - N0) * CP / (s1 - s0);
end

% Plot
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
plot(s_vec, analytical_deriv, 'b-', 'LineWidth', 2);
xlabel('s'); ylabel('dhd/ds');
title('Analytical Derivative');
grid on;

subplot(1,3,2);
plot(s_vec, numerical_deriv, 'r-', 'LineWidth', 2);
xlabel('s'); ylabel('dhd/ds');
title('Numerical Derivative');
grid on;

subplot(1,3,3);
error_vec = abs(analytical_deriv - numerical_deriv);
semilogy(s_vec, error_vec, 'k-');
xlabel('s'); ylabel('Absolute Error');
title(sprintf('Error (max = %.2e)', max(error_vec)));
grid on;

fprintf('Maximum error over full range: %.2e\n', max(error_vec));

%% Test Case 4: Check individual basis functions
fprintf('\n=== TEST CASE 4: Individual basis function derivatives ===\n');

n = 4;  % 5 control points (easier to debug)
p = 2;  % quadratic (simpler)

figure('Position', [100, 500, 1000, 300]);

s_vec = linspace(0, 1, 200);

for basis_idx = 0:n
    subplot(1, n+1, basis_idx+1);
    
    deriv_values = zeros(size(s_vec));
    for i = 1:length(s_vec)
        dN = BSpline_derivative(n, p, s_vec(i));
        deriv_values(i) = dN(basis_idx+1);
    end
    
    % Compare with numerical
    deriv_num = zeros(size(s_vec));
    for i = 1:length(s_vec)
        s = s_vec(i);
        if s < 1e-6
            s0 = s; s1 = s + 1e-6;
        elseif s > 1 - 1e-6
            s0 = s - 1e-6; s1 = s;
        else
            s0 = s - 1e-7; s1 = s + 1e-7;
        end
        N0 = BSpline(n, p, s0);
        N1 = BSpline(n, p, s1);
        deriv_num(i) = (N1(basis_idx+1) - N0(basis_idx+1)) / (s1 - s0);
    end
    
    plot(s_vec, deriv_values, 'b-', 'LineWidth', 2);
    hold on;
    plot(s_vec, deriv_num, 'r--', 'LineWidth', 1);
    xlabel('s'); 
    title(sprintf('dN_{%d,%d}/ds', basis_idx, p));
    grid on;
    
    if basis_idx == 0
        legend('Analytical', 'Numerical');
    end
end

sgtitle('Individual B-spline Basis Function Derivatives');