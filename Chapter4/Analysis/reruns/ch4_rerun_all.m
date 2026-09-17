function ch4_rerun_all(stages)
%CH4_RERUN_ALL  Every rerun the Chapter 3-4 reports need after the 2026-09-18
% changes (one actuator rating for the torque boxes, physical validity),
% in order, each stage resumable.
%
%   ch4_rerun_all()              all stages below, in order
%   ch4_rerun_all('l1')          one stage
%   ch4_rerun_all({'load','nine'})
%
% Stages, in the order they run:
%
%   robust case4 l1 load   ch4_main with that preset (25 steps, no GIFs). A
%                          preset is SKIPPED when Results/ already holds a
%                          rating-rule result for it with validity columns,
%                          so a finished preset is never rerun.
%   ch3tests               ch3_test_all (checks the Chapter 3 validity code;
%                          the collocation order test is a known failure)
%   ch3table               ch3_table_rerun: the Chapter 3 report's
%                          controller table with validity
%   long_robust long_sweep long_nine long_oos long_summary
%                          ch4_long_reruns: the robust plateau, the
%                          predictor-rate sweep, the nine long runs x
%                          variants, the six fresh load sequences, and the
%                          summary log. Each saves one file per run and skips
%                          runs already saved.
%
% Outputs: Results/ch4_result_*.mat (presets), Results/reruns/ch3 and
% Results/reruns/ch4 (the rest). All are git-ignored; to continue on another
% machine, copy them across.
%
% ONE MATLAB SESSION AT A TIME. On the author's Mac, MATLAB -batch aborts
% intermittently inside libcurl; run_reruns.sh in this folder wraps each stage
% in a retry loop that waits for MATLAB to be idle. Every stage resumes, so a
% crash only costs the run in progress.
%
% See also CH4_MAIN, CH4_LONG_RERUNS, CH3_TABLE_RERUN.

all_stages = {'robust', 'case4', 'l1', 'load', 'ch3tests', 'ch3table', ...
              'long_robust', 'long_sweep', 'long_nine', 'long_oos', 'long_summary'};
if nargin < 1 || isempty(stages), stages = all_stages; end
if ischar(stages), stages = {stages}; end

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));

for i = 1:numel(stages)
    st = stages{i};
    fprintf('\n===== ch4_rerun_all: %s =====\n', st);
    switch st
        case {'robust', 'case4', 'l1', 'load'}
            f = rating_result(root, st);
            if ~isempty(f)
                fprintf(' already done: %s\n', f);
            else
                ch4_main('presets', {st}, 'animate', false);
            end
        case 'ch3tests'
            ch3_test_all();
        case 'ch3table'
            ch3_table_rerun();
        case {'long_robust', 'long_sweep', 'long_nine', 'long_oos', 'long_summary'}
            ch4_long_reruns(st(6:end));
        otherwise
            error('ch4_rerun_all:stage', 'Unknown stage "%s".', st);
    end
    fprintf('RERUN_STAGE_DONE %s\n', st);
end
end

% ---------------------------------------------------------------------------
function f = rating_result(root, field)
% Newest Results/ch4_result_*.mat with a non-empty FIELD, run under the rating
% rule, whose rows carry the validity columns. '' if none.
f = '';
res = fullfile(root, 'Results');
F = dir(fullfile(res, 'ch4_result_*.mat'));
[~, order] = sort({F.name});
for k = fliplr(order)
    name = fullfile(res, F(k).name);
    w = whos('-file', name);
    if ~any(strcmp({w.name}, field)), continue; end
    R = load(name, field, 'p');
    rows = R.(field);
    if isempty(rows) || ~isfield(rows, 'valid_steps'), continue; end
    if ~strcmp(ch4_box_rule(R.p), 'rating'), continue; end
    f = name;
    return;
end
end
