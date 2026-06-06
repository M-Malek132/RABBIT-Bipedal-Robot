function save_results(t_all, x_all, params, impact_log, delta_log)
%SAVE_RESULTS  Save simulation results to a timestamped .mat file.
% Compatible with illustrate_results().
%
% Inputs:
%   t_all      : (N x 1) time vector
%   x_all      : (N x n_states) state trajectory (row = time)
%   params     : struct of robot parameters
%   impact_log : impact event data
%   delta_log  : (optional) constraint slack variables

    result_dir = fullfile(pwd, 'Results');
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end

    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    filename  = fullfile(result_dir, sprintf('simulation_results_%s.mat', timestamp));

    if nargin >= 5 && ~isempty(delta_log)
        save(filename, 't_all', 'x_all', 'params', 'impact_log', 'delta_log');
    else
        save(filename, 't_all', 'x_all', 'params', 'impact_log');
    end

    fprintf('Results saved to: %s\n', filename);
end
