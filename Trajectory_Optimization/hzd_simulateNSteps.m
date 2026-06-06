function [tAll, xAll, uAll] = hzd_simulateNSteps(x0, model, opt, simOpt, nSteps)
%HZD_SIMULATENSSTEPS  Simulate multiple walking steps.
%
%  After each step:
%   1. Apply impact map     : rabbit_impact_map(x_minus, params)
%   2. Apply relabelling    : rabbit_reset_map(x_plus,  params)
%
%  Adjust function names if your repo uses different names.

params = model.params;

tAll  = [];
xAll  = [];
uAll  = [];

tOffset  = 0;
xCurrent = x0;

for step = 1:nSteps

    fprintf('  Step %d / %d ... ', step, nSteps);

    try
        [t, x, u] = hzd_simulateOneStep(xCurrent, model, opt, simOpt);
        fprintf('done (%.3f s, %d pts)\n', t(end), length(t));
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        break;
    end

    tAll = [tAll;  t + tOffset];  %#ok<AGROW>
    xAll = [xAll;  x];            %#ok<AGROW>
    uAll = [uAll;  u];            %#ok<AGROW>

    % --- Impact and relabelling ---
    xMinus   = x(end,:)';
    xPlus    = rabbit_impact_map(xMinus, params);
    xCurrent = rabbit_reset_map(xPlus, params);

    tOffset  = tAll(end);
end
end
