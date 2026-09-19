function ch3_speed_up(v_targets)
%CH3_SPEED_UP  March the walking speed UP from the cold seed to 0.5 m/s.
%
%   ch3_speed_up()                      0.1232 (seed) -> 0.50 m/s
%   ch3_speed_up([0.15 0.2 0.3 0.4 0.5])
%
% WHY THIS EXISTS SEPARATELY FROM CH3_SPEED_LADDER. That driver marches DOWN
% from posture_195 and walls at 1.15 m/s: torque, stance friction and the Fz
% floor are all active there at once, and relaxing any one of them singly does
% not release it (ch3_speed_wall). Speeds below 1.15 are therefore not
% reachable by continuing that chain, and have to be approached from the other
% end -- from the hand seed, which walks at 0.1232 m/s.
%
% WHAT MAKES THIS HARDER THAN IT LOOKS. The cold seed is not a gait: its rows
% read torque +284, friction +457, impact cone +12.7. So each rung has to make
% a gait AND hold a speed, where the downward march only ever had to move a
% gait that already existed. ch3_continuation is the tool for it -- small
% speed steps keep every solve nearly feasible at its start -- and it carries
% an overshoot guard, because a step that is too big cannot be told from a
% converging one by exitflag alone.
%
% The first target is the seed's own speed, for the reason ch3_col_constraints
% gives: NEC1 must be switched on where it is already satisfied, then marched.
%
% Output: Results/reruns/speed_up/ (hist .mat per march, speed_up.log with the
% physical-row audit at every speed that converges).
%
% See also CH3_CONTINUATION, CH3_SPEED_LADDER, CH3_SPEED_WALL.

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'speed_up');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'speed_up.log');

p = ch3_params();
p.enforce_nec1  = true;
p.scale_problem = true;

z_seed = ch3_col_seed(p);
E      = ch3_col_eval(z_seed, p);
v_seed = E.L_step / E.T;

if nargin < 1 || isempty(v_targets)
    v_targets = [v_seed, 0.15, 0.175, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50];
else
    v_targets = [v_seed, v_targets(:).'];
end

logln(logf, '=== speed march UP from the cold seed | %s', datestr(now));
logln(logf, 'seed walks at %.4f m/s, N = %d nodes; targets %s', ...
      v_seed, p.N_nodes, mat2str(v_targets, 4));
report(logf, p, z_seed, 'seed', NaN);

% A cold rung has more to do than a warm one, so it gets a bigger budget than
% the downward march's rungs did (those only had to move an existing gait).
opts = struct('MaxFunctionEvaluations', 6.0e5, 'MaxIterations', 600);

[z, hist] = ch3_continuation(p, v_targets, [], z_seed, logf, opts); %#ok<ASGLU>
save(fullfile(SPD, 'march_hist.mat'), 'hist', 'p');

logln(logf, '--- audit of every speed the march reached');
for k = 1:numel(hist)
    pk = p; pk.v_des = hist(k).v_des;
    report(logf, pk, hist(k).z, sprintf('v=%.4f', hist(k).speed), NaN);
    if hist(k).verify_ok && hist(k).max_ceq < 1e-3
        z_k = hist(k).z; p_k = pk; out_k = hist(k); %#ok<NASGU>
        save(fullfile(SPD, sprintf('gait_v%04.0f.mat', 1000*hist(k).v_des)), ...
             'z_k', 'p_k', 'out_k');
    end
end
logln(logf, '=== SPEED_UP_DONE %s', datestr(now));
fprintf('SPEED_UP_DONE\n');
end

% ---------------------------------------------------------------------------
function report(logf, p, z, tag, secs)
E   = ch3_col_eval(z, p);
c   = ch3_col_constraints(z, p);
chk = ch3_col_check_limits(z, p);
V   = ch3_col_verify(z, p, false);
pk  = max([abs(E.u(:)); abs(E.um(:))]);
lam = [E.lam, E.lamm];
mu_step = max(abs(lam(1,:)) ./ lam(2,:));
mu_land = abs(E.impulse(1)) / E.impulse(2);
logln(logf, ['%-11s v %.4f m/s | T %.4f | L %.4f | peak|u| %8.2f Nm (%.0f%% of box) | ' ...
             'mu step %.3f land %.4f | minFz %.1f N | limits ok %d (max|c| %.2e) | ' ...
             'verify %d (%.1e) | %.0f s'], ...
      tag, E.L_step / E.T, E.T, E.L_step, pk, 100 * pk / p.limits.u_max, ...
      mu_step, mu_land, min(lam(2,:)), chk.ok, chk.max_c, V.ok, V.max_dev, secs);
names = {4,'torque';5,'stance friction';6,'Fz floor';13,'impact cone';9,'swing clearance';14,'HZD exists';15,'HZD stable'};
bad = '';
for i = 1:size(names,1)
    k = names{i,1};
    if c(k) > 1e-6, bad = [bad sprintf(' %s(+%.1e)', names{i,2}, c(k))]; end %#ok<AGROW>
end
if isempty(bad), bad = ' none'; end
logln(logf, '      violated:%s', bad);
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
