function lim = ch5_robust_ulim(runs)
%CH5_ROBUST_ULIM  A shared control-axis limit that one spike cannot destroy.
%
%   lim = ch5_robust_ulim(runs)
%
% The control axis is shared across columns on purpose -- cross-column
% comparison is the claim both Fig. 5.3 and Fig. 5.4 are making. But sharing an
% axis makes it hostile to outliers, and these runs have exactly one kind.
%
% The min-norm CLF-QP applies ||mu|| = psi/||LgV||, which is unbounded wherever
% LgV passes near zero. On the pendulum baseline that happens once, at
% t = 0.224 s, where ||LgV|| dips 140x below its median and the torque touches
% 435 Nm for a single sample against a 99th percentile of 27. Scaling the
% shared axis to 435 would compress every real feature in all three columns
% into a flat line to accommodate 1 sample in 30001.
%
% So the limit comes from the 99.5th percentile across all runs, and any panel
% whose data leaves the axis is ANNOTATED WITH ITS TRUE PEAK by
% ch5_note_clipping. Clipping silently would be misleading; clipping with the
% number printed on the panel is just a readable axis.
%
% See also CH5_NOTE_CLIPPING, CH5_PLOT_PENDULUM, CH5_PLOT_SPRINGMASS.

all_u = [];
for i = 1:numel(runs)
    all_u = [all_u, max(abs(runs(i).sim.u), [], 1)]; %#ok<AGROW>
end
all_u = all_u(isfinite(all_u));

if isempty(all_u)
    lim = 1;
    return;
end

lim = 1.25 * max(ch5_prctile(all_u, 99.5), eps);

% If nothing is actually an outlier, show everything -- do not crop a clean
% signal just because a percentile said so.
if max(all_u) <= 2 * lim
    lim = 1.1 * max(all_u);
end

end
