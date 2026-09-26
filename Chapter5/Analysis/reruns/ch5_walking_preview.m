function A = ch5_walking_preview()
%CH5_WALKING_PREVIEW  Corollary 5.2 at the walking robot's post-impact state.
%
%   A = ch5_walking_preview()
%
% Chapter 5 validates the ECBF on two plants that start at rest, where the
% admissibility of the poles (Corollary 5.2) holds with full margin, and names
% the dependence on x0 as its main limitation. On a walking robot that
% limitation is not a footnote: every step starts from the output of the
% previous impact, so the condition has to hold again at every step. This
% makes the hand-off concrete. It takes the reference gait's fixed point --
% posture_195, the gait Chapters 3, 4 and 6 run on -- and its post-impact
% state x+ (the start of the periodic step), and evaluates Chapter 6's
% footstep barriers there with Chapter 6's own defaults (ch6_params: the
% stepping-stone problem, the nominal stone):
%
%   g(x+)          the barrier itself            (C_0: must be >= 0)
%   gdot(x+)       its rate at the start of the step
%   gamma_b        the first pole Chapter 6 uses for it
%   gamma_b,min    = -gdot / g, the smallest first pole Corollary 5.2 allows
%   h_CBF(x+)      = gamma_b g + gdot (C_1: must be >= 0), and its margin
%
% A barrier whose gamma_b,min is close to the pole in use is one the next
% disturbed impact can make inadmissible. Chapter 6 rechecks this at every
% step (ch6_step records ch6_admissible for each one); this is the number at
% the undisturbed orbit, i.e. the margin those checks start from.
%
% Output: struct A, and Results/reruns/ch5_walking_preview/preview.log.
% Runtime: seconds.
%
% See also CH6_ADMISSIBLE, CH6_BARRIER, CH6_LOAD_GAIT, CH5_ECBF_ADMISSIBLE.

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
OUT  = fullfile(ROOT, 'Results', 'reruns', 'ch5_walking_preview');
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, 'preview.log');
if exist(logf, 'file'), delete(logf); end

[x0, alpha, p6, meta] = ch6_load_gait();
a = ch6_admissible(x0, alpha, p6, 0);
[B, ~] = ch6_barrier(x0, [], p6, 0);

logln(logf, '=== ch5_walking_preview | %s | %s', meta.file, datestr(now)); %#ok<TNOW1,DATST>
logln(logf, 'post-impact state of the fixed point: step %.4f m in %.4f s (%.3f m/s)', ...
      meta.L_step, meta.T, meta.v_avg);
A = struct('label', {}, 'g', {}, 'gdot', {}, 'gamma_b', {}, 'gamma_b_min', {}, ...
           'h_cbf', {}, 'ok', {});
for i = 1:numel(B)
    gb = poles_for(p6, B(i).label);
    g  = B(i).g;  gd = B(i).gdot;
    gmin = NaN;
    if g > 0, gmin = max(0, -gd / g); end
    hc = gb * g + gd;
    A(end+1) = struct('label', B(i).label, 'g', g, 'gdot', gd, 'gamma_b', gb, ...
                      'gamma_b_min', gmin, 'h_cbf', hc, 'ok', g >= 0 && hc >= 0); %#ok<AGROW>
    logln(logf, ['%-16s g %+.4f | gdot %+.4f | gamma_b %.1f, admissible down to %.2f ' ...
                 '(margin x%.1f) | h_CBF %+.4f | %s'], B(i).label, g, gd, gb, gmin, ...
          gb / max(gmin, eps), hc, ternary(A(end).ok, 'admissible', 'NOT admissible'));
end
logln(logf, 'ch6_admissible: ok %d%s', a.ok, ternary(a.ok, '', [' -- ' a.why]));
save(fullfile(OUT, 'preview.mat'), 'A', 'a');
logln(logf, '=== WALKING_PREVIEW_DONE');
end

% ---------------------------------------------------------------------------
function gb = poles_for(p, label)
% The first pole ch6_cbf_row uses for this barrier: its own pair from
% p.cbf.poles_by_label when listed there, else p.cbf.gamma_b.
gb = p.cbf.gamma_b;
if isfield(p.cbf, 'poles_by_label') && ~isempty(p.cbf.poles_by_label)
    L = p.cbf.poles_by_label;
    for r = 1:size(L, 1)
        if strcmp(L{r, 1}, label), gb = L{r, 2}(1); end
    end
end
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
