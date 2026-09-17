function f = ch4_figure(name, width, height_cm)
%CH4_FIGURE  A Chapter-4 figure, drawn at the size the report prints it.
%
%   f = ch4_figure(name, width, height_cm)
%
% width is 'column' (8.4 cm, one column of the report's two-column A4 page) or
% 'page' (17.4 cm, across both columns).
%
% DRAWN AT PRINT SIZE, TEXT IN POINTS. The report scales every figure to the
% width it is set at, so whether its text can be read depends on the text's
% size relative to the figure, not on the resolution. These figures used to be
% drawn in a 22 cm window at MATLAB's 10 pt default; set in an 8.4 cm column
% that is 3-4 pt, and no title or axis in the report could be read. Drawn at
% 8.4 cm with 9 pt text, the text prints at 9 pt.
%
% Save with exportgraphics, which crops to the content, rather than saveas,
% which keeps the window's margins and renders at 150 dpi.
%
% See also CH4_PLOT_UNCERTAINTY, CH4_PLOT_LOAD, CH4_DOC_FIGURES.

widths = struct('column', 8.4, 'page', 17.4);

% Legends: no box, whose white background, on a legend set inside a panel,
% hid the peaks of the curves under it; and short line samples, since at the
% default length a three-entry legend was wider than a column.
f = figure('Name', name, 'Color', 'w', 'Units', 'centimeters', ...
           'Position', [1 1 widths.(width) height_cm], ...
           'DefaultAxesFontSize', 9, ...               % ticks; legends get 90%
           'DefaultAxesLabelFontSizeMultiplier', 1.1, ...
           'DefaultAxesTitleFontSizeMultiplier', 1.1, ...
           'DefaultTextFontSize', 9, ...
           'DefaultAxesLineWidth', 0.5, ...
           'DefaultLineLineWidth', 0.8, ...
           'DefaultLegendBox', 'off', ...
           'DefaultLegendItemTokenSize', [15 8]);
end
