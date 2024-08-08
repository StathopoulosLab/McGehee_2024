function [stat1,p1, stat2, p2] = calc_start_var(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [stat1,p1] = grouped_plot(data);
    [stat2,p2] = ungrouped_plot(data);
end

function [stat,p] = grouped_plot(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    % Initialize looped variable
    n = zeros(size(data, 2), 1);
    c = zeros(size(data, 2), 1);
    g = strings(size(data, 2), 1);

    unique_conditions = {'dark', '10', '20'}; %unique({data.condition}, 'stable');
    
    % For each data set
    for i = 1:size(data, 2)
        % Calculate number of detected spots for the first time point
        n(i) = size(data(i).A{1,1},1);
        temp = strsplit(data(i).condition, ' ');
        c(i) = find(strcmp(temp{2},unique_conditions));
        g(i) = temp{1};
    end

    [stat, p] = plot_val(n, g, c, [2,1], [8,1,4], {'{\itdl-LEXY}', '{\itdl-BLID}'}, true);
end

function [stat,p] = ungrouped_plot(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    % Initialize looped variable
    n = zeros(size(data, 2), 1);
    c = zeros(size(data, 2), 1);
    g = strings(size(data, 2), 1);

    unique_conditions = unique({data.condition}, 'stable');
    
    % For each data set
    for i = 1:size(data, 2)
        % Calculate number of detected spots for the first time point
        n(i) = size(data(i).A{1,1},1);
        c(i) = find(strcmp(data(i).condition,unique_conditions));
        g(i) = data(i).condition;
    end

    [stat, p] = plot_val(n, g, c, [5,6,7,1,2,3], [8,1,4,8,1,4], {'{\itdl-LEXY}', '{\itdl-BLID}'}, false);
end

function [stat, p] = plot_val(y, g, c_ind, x_order, c_order, x_tick_labels, grouped)
%PLOT_MEAN_INTENSITY Plot the individual intensity, mean, and error bars
% 
%   Inputs
%       avg: the avg structure from calc_means
%       field: the field to be plotted
% 
%   Outputs
%       None
% 
%   Overview
%       A plot is made of the individual data, the mean, and the error bars
%       for the data specified by field
            
    % Array of colors
    colors = [     0, 0.4470, 0.7410;   %blue
              0.8500, 0.3250, 0.0980;   %red
              0.9290, 0.6940, 0.1250;   %yellow/gold
              0.4940, 0.1840, 0.5560;   %purple
              0.4660, 0.6740, 0.1880;   %green
              0.3010, 0.7450, 0.9330;   %cyan
              0.6350, 0.0780, 0.1840;   %dark red
                   0,      0,      0];  %black
    
    % Save condition name for legend
    legend_names = unique(g, 'stable');
%     x = strcmp(g, legend_names{2}) + 1;
    
    x = zeros(size(y, 1), 1);
    mean_y = zeros(size(legend_names, 1), 1);
    err_y = zeros(size(legend_names, 1), 1);
    c = zeros(size(y, 1), 3);

    for i = 1:size(legend_names, 1)
        x(strcmp(g,legend_names{i})) = x_order(i);
        mean_y(i) = mean(y(strcmp(g,legend_names{i})), 1);
        err_y(i) = calc_error(y(strcmp(g,legend_names{i})), 'SD', 1);

        for j = 1: size(y,1)
            c(j,:) = colors(c_order(c_ind(j)),:);
            %c(strcmp(g,legend_names{i}),:) = repmat(colors(i,:), [sum(strcmp(g,legend_names{i}), 1), 1]);
        end
    end
    
    [stat, p] = statistical_analysis(y, g);

    % Make figure
    figure;

    % Set formating for labels
    set(gca, 'TickLabelInterpreter', 'tex');

    % Allow plotting without overwriting
    hold on;

    % Plot scatter plot of individual data
    scatter(x, y, 500, c, 'filled',...
            'jitter', 'on', 'jitteramount', 0.2,...
            'MarkerFaceAlpha', .85,'MarkerEdgeAlpha', .85);
    scatter(nan,nan, 500, colors(c_order(2),:), 'filled');
    scatter(nan,nan, 500, colors(c_order(3),:), 'filled');
    
    % Plot the mean with the error bars and set properties
    h = errorbar(x_order, mean_y, err_y, '.', 'Color', [0.5, 0.5, 0.5]);
    set(h, 'linewidth', 4, 'markersize', 50);
    % Turn off plotting without overwriting
    hold off;

    % Set properties of axis
        set(gca, 'ylim', [0, 400],...
                 'ytick', 0:100:400,...
                 'xlim', [0, max(x_order)+1],...
                 'XTickLabelRotation', 0,...
                 'fontsize', 30,...
                 'fontname', 'arial');
    
    if grouped
        set(gca, 'xtick', 1:1:size(mean_y,1));%((size(avg, 2)+1)/2):size(avg, 2):size(y, 2),...)
    else
        set(gca, 'xtick', 2:4:size(mean_y,1));
    end

    % Label x ticks
%         xticklabels(x_tick_labels(nc));
    xticklabels(x_tick_labels);

    % Label y
%         ylabel(y_labels{j}, 'FontSize', 18);
    % Make a legend
    special_legend_names = {'Dark', 'Before light 10 min', 'Before light 20 min'};
    [~, hobj, ~, ~] = legend(special_legend_names{:}, 'Location', 'northwest', 'Fontsize', 20, 'fontname', 'arial');
    legend('boxoff');
    M = findobj(hobj,'type','patch');
    set(M,'MarkerSize',sqrt(500));
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

function [stat, p] = statistical_analysis(data, condition)
%STATISTICAL_ANALYSIS Perform ANOVA to compare the means between conditions
% 
%   Input
%       avg: the structure returned from calc_means
% 
%   Output
%       stat: structure containing outputs from anova1 and multcompare
%       p: table of p-values for pairwise comparisons
% 
%   Overview
%       This function performs statistical analysis on the data.
%       Specifically, it performs one way ANOVA using anova1 and multiple
%       comparisons using Tukey's HSD using multcompare. It returns the
%       outputs from anova1 and multcompare in the structure stat and
%       a table of p-values, p, for pairwise comparison between conditions.
    
%     condition = cat(1, avg(1,:).condition);
%     data = cat(1, avg(1,:).(field));
    
    % Initialize a structure for storing the results of the statistical
    % analysis
    stat = struct('p', [],...
                  'tbl', [],...
                  'stats', [],...
                  'p_indiv', [],...
                  'means', [],... 
                  'names', []);
    
    % Perform ANOVA on the intensity data grouped by condition
    [stat.p, stat.tbl, stat.stats, stat.p_indiv, stat.means,...
        stat.names] = stat_test(data, condition);
    
    % Initialize variables for making comparison tables
    p = cell(size(stat.names, 1), size(stat.names, 1));

    % Save time point i in table of mutiple comparisons
    p{1, 1} = 'p-values';
    
    % Make row names of conditions for comparison
    p(2:end, 1) = stat.names(1:(end-1));
    
    % Make column names of conditions for comparison
    p(1, 2:end) = stat.names(2:end);
            
    % For each comparison
    for j = 1:size(stat.p_indiv, 1)
        % save the p-value in the p-value table
        p{stat.p_indiv(j,1) + 1,...
            stat.p_indiv(j,2)} = stat.p_indiv(j,6);
    end
end

function [p, tbl, stats, p_indiv, means, names] = stat_test(data,...
    group)
%STAT_TEST Perform ANOVA to compare the means between conditions
% 
%   Input
%       data: data that anova will be performed on
%       group: identifier for data to correctly group it
% 
%   Output
%       p: p-value from the anova
%       tbl: a table returned from anova
%       stats: statistics for mutiple comparison tests
%       p_indiv: pairwise p-values from mutiple comparisons
%       means: estimated means
%       names: names of groups
% 
%   Overview
%       This function performs statistical analysis on the data. 
%       Specifically, it performs one way ANOVA using anova1 and multiple
%       comparisons using Tukey's HSD using multcompare. It returns the
%       outputs from anova1 and multcompare in the structure stat.

    % Perform ANOVA on the data grouped by condition in group
    [p, tbl, stats] = anova1(data, group, 'off');
    
    % Perform pairwise comparisons of data between conditions
    % using Tukey's HSD
    [p_indiv, means, ~, names] = multcompare(stats, 'display', 'off');
end