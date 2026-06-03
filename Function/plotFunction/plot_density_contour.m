function plot_density_contour(x, y, gridpts, xg, yg, color, xlimv, ylimv)

    keep = ~isnan(x) & ~isnan(y);
    data = [x(keep) y(keep)];

    f = ksdensity(data, gridpts);
    F = reshape(f, size(xg));
    % 
    % contourf(xg, yg, F, 8, 'LineColor','none');
    % hold on
    % contour(xg, yg, F, 6, 'LineColor', color, 'LineWidth',1.2);
    % 
    % plot([0 600],[0 600],'--','Color',[0.6 0.6 0.6],'LineWidth',1)
    % 
    % axis image
    % % xlim([0 600]); ylim([0 600])
    % colormap(gca, flipud(gray))

    n = 256;
    cmap = [linspace(1,color(1),n)', ...
            linspace(1,color(2),n)', ...
            linspace(1,color(3),n)'];

    contourf(xg, yg, F, 10, 'LineColor','none');
    hold on
    contour(xg, yg, F, 6, 'LineColor', color, 'LineWidth',1.3);

    colormap(gca, cmap)

    plot([0 ylimv],[0 ylimv],'--','Color',[0.6 0.6 0.6],'LineWidth',1)

    axis image
    xlim([0 xlimv]); ylim([0 ylimv])
    box off; set(gca,'TickDir','out','LineWidth',1.2,'FontSize',11)
end