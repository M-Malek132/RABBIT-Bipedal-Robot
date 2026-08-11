function pass = ch6_test_library()
%CH6_TEST_LIBRARY  The gait library map of (6.22)-(6.24).
%
% This suite runs on a SYNTHETIC library -- coefficients that are a known
% function of step length -- rather than on solved gaits. That is deliberate:
% the claim being checked is about the interpolation, and with a synthetic
% library the exact answer is known, so a wrong bracket or a reversed zeta is
% caught by a number rather than by a plot looking odd. Whether the solved
% gaits are good gaits is ch6_lib_build's business and ch3_col_verify's.
%
% Checks:
%   1. at a grid point, alpha is that gait exactly (zeta = 0 or 1)
%   2. between grid points, (6.23) reproduces the linear ground truth
%   3. outside the grid, it EXTRAPOLATES linearly and says so
%   4. the bracket is the straddling pair, including at the last interval
%   5. a real library on disk, if present, loads and is ascending
%
% See also CH6_LIB_ALPHA, CH6_LIB_BUILD, CH6_LIB_LOAD.

fprintf('\n=== ch6_test_library ===\n');
pass = true;

ny = 4;  nc = 6;
L  = [0.10 0.24 0.40 0.56 0.72];         % the thesis's own grid shape
m  = numel(L);

% Ground truth: alpha(L) = A0 + A1*L, exactly linear, so interpolation must be
% exact everywhere including outside the grid.
A0 = randn(ny, nc);
A1 = randn(ny, nc);
lib = struct('L', L, 'alpha', zeros(ny, nc, m));
for i = 1:m
    lib.alpha(:,:,i) = A0 + A1 * L(i);
end

%% ---------------------------------------------------- 1. exact at grid points
err = 0;
for i = 1:m
    [a, info] = ch6_lib_alpha(lib, L(i));
    err = max(err, norm(a - lib.alpha(:,:,i), 'fro'));
    err = max(err, abs(info.zeta - (i == m)));   % last point: zeta = 1
end
ok = err < 1e-14;
fprintf('  [%s] %-42s max err %.1e\n', tf(ok), ...
        'alpha(L_i) is gait i exactly', err);
pass = pass && ok;

%% ------------------------------------------ 2. exact between grid points
err = 0;
for Ld = linspace(L(1), L(end), 97)
    a = ch6_lib_alpha(lib, Ld);
    err = max(err, norm(a - (A0 + A1*Ld), 'fro'));
end
ok = err < 1e-13;
fprintf('  [%s] %-42s max err %.1e\n', tf(ok), ...
        '(6.23) reproduces a linear family', err);
pass = pass && ok;

%% ---------------------------------------------- 3. extrapolation, reported
err = 0;  flagged = true;
for Ld = [0.02 0.05 0.90 1.10]
    [a, info] = ch6_lib_alpha(lib, Ld);
    err = max(err, norm(a - (A0 + A1*Ld), 'fro'));
    flagged = flagged && info.extrapolated;
end
[~, in] = ch6_lib_alpha(lib, 0.30);
ok = err < 1e-13 && flagged && ~in.extrapolated;
fprintf('  [%s] %-42s err %.1e, flagged %d\n', tf(ok), ...
        'linear extrapolation past the ends', err, flagged);
pass = pass && ok;

% and the amount of extrapolation, in units of the end interval (Remark 6.6)
[~, i9] = ch6_lib_alpha(lib, 0.72 + 0.5*(0.72-0.56));
ok = abs(i9.extrap_frac - 0.5) < 1e-12;
fprintf('  [%s] %-42s extrap_frac = %.4f (expect 0.5)\n', tf(ok), ...
        'extrapolation distance is reported', i9.extrap_frac);
pass = pass && ok;

%% ---------------------------------------------------- 4. the bracket straddles
bad = 0;
for Ld = linspace(L(1), L(end), 200)
    [~, info] = ch6_lib_alpha(lib, Ld);
    if ~(L(info.i) <= Ld + 1e-12 && Ld <= L(info.i+1) + 1e-12), bad = bad + 1; end
end
ok = bad == 0;
fprintf('  [%s] %-42s %d/200 wrong brackets\n', tf(ok), ...
        'the bracket straddles the request', bad);
pass = pass && ok;

%% -------------------------------------------------- 5. a real library, if any
p = ch6_params();
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
f = fullfile(root, 'Results', 'ch6_gait_library.mat');
if exist(f, 'file')
    try
        rl = ch6_lib_load(p);
        ok = numel(rl.L) >= 1 && all(diff(rl.L) > 0);
        fprintf('  [%s] %-42s %d gaits, L in [%.3f, %.3f]\n', tf(ok), ...
                'the library on disk loads and is sorted', ...
                numel(rl.L), rl.L(1), rl.L(end));
    catch ME
        ok = false;
        fprintf('  [%s] %-42s %s\n', tf(ok), ...
                'the library on disk loads', ME.message);
    end
    pass = pass && ok;
else
    fprintf('  [skip] %-40s no Results/ch6_gait_library.mat yet\n', ...
            'library on disk');
end

fprintf('  --> %s\n', upper(tf(pass)));

end

function s = tf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
