function plot_scatter_panel(x, y, color, ms, alpha, xlimv, ylimv)

    keep = ~isnan(x) & ~isnan(y);
    x = x(keep); y = y(keep);

    scatter(x, y, ms, color, 'filled', ...
        'MarkerFaceAlpha', alpha, ...
        'MarkerEdgeAlpha', alpha)
    hold on

    plot([0 ylimv], [0 ylimv], '--', 'Color',[0.6 0.6 0.6], 'LineWidth',1)

    axis image
    xlim([0 xlimv]); ylim([0 ylimv])

end
