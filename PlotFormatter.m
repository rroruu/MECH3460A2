function PlotFormatter()

    reset(groot);

    set(groot, 'defaultFigureUnits', 'centimeters');
    set(groot, 'defaultFigurePosition', [2, 2, 8.89, 5.5]);
    set(groot, 'defaultFigureColor', 'w');

    set(groot, 'defaultLineLineWidth', 1.2);
    set(groot, 'defaultLineColor', [0 0.4470 0.7410]); 
    set(groot, 'defaultAxesFontName', 'Times New Roman');
    set(groot, 'defaultAxesFontSize', 8); 
    set(groot, 'defaultAxesBox', 'off');
    set(groot, 'defaultAxesXGrid', 'on');
    set(groot, 'defaultAxesYGrid', 'on');
    set(groot, 'defaultAxesGridAlpha', 0.15);

    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    set(groot, 'defaultLegendLocation', 'northoutside');
    set(groot, 'defaultLegendOrientation', 'horizontal');
    set(groot, 'defaultLegendBox', 'off');
    set(groot, 'defaultLegendFontSize', 7);

    set(groot, 'defaultFigureCreateFcn', @(fig, ~) addlistener(fig, 'SizeChanged', ...
        @(src, ~) arrayfun(@(ax) set(ax, 'LooseInset', get(ax, 'TightInset')), ...
        findobj(src, 'Type', 'axes'))));

    disp('IEEE plot defaults successfully set for this session.');
end