function f = ch5_figure(name, width, height_cm)
%CH5_FIGURE  A Chapter-5 figure, drawn at the size the report prints it.
%
%   f = ch5_figure(name, width, height_cm)
%
% width is 'column' (8.4 cm, one column of the report's two-column A4 page) or
% 'page' (17.4 cm, across both columns).
%
% DRAWN AT PRINT SIZE, TEXT IN POINTS. The report scales every figure to the
% width it is set at, so whether its text can be read depends on the text's
% size relative to the figure, not on the resolution. These figures used to be
% drawn in a 26-35 cm window at MATLAB's 10 pt default and set across a 17.4 cm
% page, which printed their titles and axes at about 5 pt.
%
% The sizes and fonts are ch4_figure's, so the two reports' figures match. Save
% with exportgraphics at 600 dpi.
%
% See also CH4_FIGURE, CH5_PLOT_SPRINGMASS, CH5_PLOT_PENDULUM, CH5_DOC_FIGURES.

widths = struct('column', 8.4, 'page', 17.4);

% Legends: no box, whose white background hides the curves under a legend set
% inside a panel; and short line samples, or a legend outgrows its panel.
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
