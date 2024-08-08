function data = plot_FRAP_quant(data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Array of colors
    colors = [0.4660, 0.6740, 0.1880;   %green
              0.3010, 0.7450, 0.9330;   %blue
              0.9290, 0.6940, 0.1250;   %yellow/gold
              0.8500, 0.3250, 0.0980;   %red
              0.4940, 0.1840, 0.5560];   %purple];   %cyan

    figure;
    hold on;
    
    for i = 1:size(data, 2)
%         norm = data(i).avg_I(20);
        norm = 1;
        a0 = [mean(data(i).avg_I((end-10):end))/norm, min(data(i).avg_I)/norm, 1];
        data(i).params = fit_recovery_curve(a0, data(i).time(data(i).t_bleach:end)-data(i).time(data(i).t_bleach),...
            data(i).avg_I(data(i).t_bleach:end)/norm);
    
%         h = plot(data(i).time, data(i).avg_I, 'Color', colors(i,:));
%         set(h, 'linewidth', 1.5);

        scatter(data(i).time, data(i).avg_I/norm, 40, colors(i,:), 'filled',...
            'MarkerFaceAlpha', .5,'MarkerEdgeAlpha', .5);

        N = data(i).params(1) - data(i).params(2) * exp(-(data(i).time(data(i).t_bleach:end)-data(i).time(data(i).t_bleach)) * data(i).params(3));
        h1 = plot(data(i).time(data(i).t_bleach:end), N, 'Color', 1 * colors(i,1:3));
        set(h1, 'linewidth', 4);
    end

    hold off;

    set(gca, 'xlim', [0,23]);%[0,10]);
    set(gca, 'ylim', [0,12000]);%[0,5000]);[0,12000])[0,1.1]
    set(gca, 'Fontsize', 20);
    set(gca, 'fontname', 'arial');
end

function af = fit_recovery_curve(a0, t, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    f = @(a) a(1) - a(2) * exp(-a(3) * t) - y;

    options = optimset('Display','off');
    ub = [Inf, Inf, Inf];
    lb = [0, 0, 0];
    af = lsqnonlin(f, a0, lb, ub, options);
end       