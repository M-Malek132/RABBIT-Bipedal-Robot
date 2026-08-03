function ch5_note_clipping(ax, sim, lim, unit)
%CH5_NOTE_CLIPPING  Print the true peak on any panel whose data left the axis.
%
%   ch5_note_clipping(ax, sim, lim, unit)
%
% A shared axis chosen robustly (ch5_robust_ulim) will crop the occasional
% min-norm spike. Cropping is the right call for readability and the wrong call
% for honesty unless the reader is told, so this stamps the actual peak and how
% many samples reached it onto the panel itself.
%
% Nothing is drawn when the data fits.
%
% See also CH5_ROBUST_ULIM, CH5_PLOT_PENDULUM, CH5_PLOT_SPRINGMASS.

au = max(abs(sim.u), [], 1);
pk = max(au);

if ~isfinite(pk) || pk <= lim
    return;
end

n_out  = sum(au > lim);
[~, k] = max(au);

% The note goes in the panel TITLE, not in free-floating text. Normalized text
% coordinates are resolved against whatever axes are current when the figure is
% finally drawn, and with subplot panels that get resized afterwards -- row 3
% here calls axis equal -- the label can end up rendered against a different
% panel entirely. A title is anchored to the axes object and cannot drift.
title(ax, sprintf('peak %.0f %s at t = %.2f s (%d sample%s off-axis)', ...
                  pk, unit, sim.t(k), n_out, plural(n_out)), ...
      'FontSize', 8, 'FontWeight', 'normal', 'Color', [0.6 0 0]);

end

% ---------------------------------------------------------------------------
function s = plural(n)
if n == 1, s = ''; else, s = 's'; end
end
