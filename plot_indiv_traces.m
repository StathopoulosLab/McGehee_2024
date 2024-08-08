function data = plot_indiv_traces(data)
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
        data(i).params_on = zeros(size(data(i).I,1),3);
        data(i).params_off = zeros(size(data(i).I,1),3);
        
        for j = 1:size(data(i).I,1)
%             norm = mean(data(i).I(j,1:10), 'omitnan');
            a0_on = [min(data(i).I(j,:)), min(data(i).I(j,:)) - mean(data(i).I(j,1:10), 'omitnan'), 1];
            t_on = data(i).time(data(i).t_blue_on:data(i).t_blue_off)-data(i).time(data(i).t_blue_on);
            I_on = data(i).I(j,data(i).t_blue_on:data(i).t_blue_off);

            data(i).params_on(j,:) = fit_recovery_curve(a0_on, t_on(~isnan(I_on)),...
                I_on(~isnan(I_on))');

            a0_off = [mean(data(i).I(j,(end-10):end), 'omitnan'), mean(data(i).I(j,(end-10):end), 'omitnan') - min(data(i).I(j,:)), 1];
            t_off = data(i).time(data(i).t_blue_off:end)-data(i).time(data(i).t_blue_off);
            I_off = data(i).I(j,data(i).t_blue_off:end);

            data(i).params_off(j,:) = fit_recovery_curve(a0_off, t_off(~isnan(I_off)),...
                I_off(~isnan(I_off))');

            scatter(data(i).time, data(i).I(j,:), 40, colors(i,:), 'filled',...
                'MarkerFaceAlpha', .05,'MarkerEdgeAlpha', .05);
        end

%         data(i).mean_params_on = cat(1, mean(data(i).params_on), calc_error(data(i).params_on, 'SD', 1));
%         N1 = data(i).mean_params_on(1,1) - data(i).mean_params_on(1,2)...
%             * exp(-(data(i).time(data(i).t_blue_on:data(i).t_blue_off)-data(i).time(data(i).t_blue_on)) * data(i).mean_params_on(1,3));
%         h1 = plot(data(i).time(data(i).t_blue_on:data(i).t_blue_off), N1, 'Color', colors(i,:));
%         set(h1, 'linewidth', 4);

        data(i).mean_params_off = cat(1, mean(data(i).params_off), calc_error(data(i).params_off, 'SD', 1));
        N2 = data(i).mean_params_off(1,1) - data(i).mean_params_off(1,2)...
            * exp(-(data(i).time(data(i).t_blue_off:end)-data(i).time(data(i).t_blue_off)) * data(i).mean_params_off(1,3));
        h2 = plot(data(i).time(data(i).t_blue_off:end), N2, 'Color', colors(i,:));
        set(h2, 'linewidth', 4);
    end

    hold off;

    set(gca, 'xlim', [0, 31]);
    set(gca, 'ylim', [2000, 7000]);%[2000, 7000][0.4, 1.1]
    set(gca, 'Fontsize', 20);
    set(gca, 'fontname', 'arial');
end

function af = fit_recovery_curve(a0, t, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    f = @(a) a(1) - a(2) * exp(-a(3) * t) - y;

    options = optimset('Display','off');
    ub = [Inf, Inf, Inf];
    lb = [0, -Inf, -Inf];
    af = lsqnonlin(f, a0, lb, ub, options);
end

function err_d = calc_error(d, select_error, dim)
%CALC_EROR Calculates the error for data d
% 
%   Input
%       d: data points 
%       select_error: a string, either SEM, CI, or SD
% 
%   Output
%       err_d: the error for the data
% 
%   Overview
%       This function calculates the error for determining error bars. It
%       takes data d and the choice for calculating the error, either
%       standard error of the mean (SEM), 95% confidence intervals (CI), or
%       standard deviation (SD)
    
    % standard deviation for data in d
    STD_d = std(d, [], dim, 'omitnan');

    % standard error of the mean for data in d
    SEM_d = STD_d ./ sqrt(sum(~isnan(d), dim));

    % confidence interval for data in d
    ts_d = tinv(0.975, sum(~isnan(d), dim) - 1);
    CI_d = ts_d .* SEM_d;

    % If user inputed SD
    if ~isempty(select_error) && isequal(select_error, 'SD')
        % Error is standard deviation
        err_d = STD_d;
    % Else if user inputed CI
    elseif ~isempty(select_error) && isequal(select_error, 'CI')
        % Error is confidence intervals
        err_d = CI_d;
    % Else if user inputed SEM or anything else
    else
        % Error is standard error of the mean
        err_d = SEM_d;
    end
end