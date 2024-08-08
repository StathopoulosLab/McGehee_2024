function [stat,p,avg] = plot_gene_props(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    % Find the unique conditions in data
%     conditions = unique({data{i}.condition}', 'stable');
    
    % Calulate the mean
    avg = calc_means(data);
    [stat, p] = plot_width(avg);
    

%     [stat, p, p_ttest] = statistical_analysis(avg, field, false);
    
%     if ~norm_true
%         % Plot the mean for the given field
%         plot_mean_intensity(avg{i}, field, x_order_label{i}, c_order{i},...
%             norm_true, control_cond{i});
%     end
% 
%     if norm_true
%         % Plot the mean for the given field
%         plot_mean_intensity(avg, field, x_order_label{:}, c_order{:},...
%             norm_true, control_cond);
%     end
end

function tidy_tbl = calc_means(data)
%CALC_MEANS Calculate the mean intensity
% 
%   Inputs
%       data: data structure from main function from quantify_in_situ
%       condition: the condition of treatment
% 
%   Outputs
%       avg: a structure containg the average data, the individual data,
%           and the error bar data
% 
%   Overview
%       Calculates the mean and error for data

    % Array of colors
    colors = [     0, 0.4470, 0.7410;   % blue
              0.8500, 0.3250, 0.0980;   % red
              0.4660, 0.6740, 0.1880;   % green
              0.8500, 0.3250, 0.0980;   % red
              0.4940, 0.1840, 0.5560;   % purple
              0.9290, 0.6940, 0.1250;   % yellow/gold
              0.3010, 0.7450, 0.9330;   % cyan
              0.6350, 0.0780, 0.1840;   % dark red
                 0.5,    0.5,    0.5;   % gray
              0.5882, 0.2941,     0];   % brown
            
    init_zeros = zeros(sum(cat(2,data.n_domains)), 1);
    init_zeros2 = zeros(sum(cat(2,data.n_domains)), 3);
    init_string = strings(sum(cat(2,data.n_domains)), 1);
    init_tidy = {init_string, init_zeros, init_string, init_string,init_zeros,init_zeros2};

    % Make a table using the cell array, so that there are t rows
    % and 2 column arrays in each table variable
    tidy_tbl = table(init_tidy{:},...
                    'VariableNames', {'gene', 'gene_width', 'condition', 'stage', 'x', 'color'});

    ind = 1;
               
    % For each entered condition
    for i = 1:size(data, 2)
        for j = 1:size(data(i).gene_length,1)
            for k = 1:size(data(i).gene_length{j},1)
                tidy_tbl.gene(ind) = data(i).gene_id{j}{k};
                tidy_tbl.gene_width(ind) = data(i).gene_length{j}(k);
                tidy_tbl.condition(ind) = data(i).condition;
                tidy_tbl.stage(ind) = data(i).stage;

                unique_cond = unique(tidy_tbl.condition, 'stable');
                x = 1:size(unique_cond,1);

                tidy_tbl.x(ind) = x(strcmp(tidy_tbl.condition(ind), unique_cond));

                unique_gene = unique(tidy_tbl.gene, 'stable');
                c = 1:size(unique_gene,1);
                tidy_tbl.color(ind,:) = colors(c(strcmp(tidy_tbl.gene(ind), unique_gene)),:);

                if strcmp(tidy_tbl.gene(ind), 'zen') && strcmp(tidy_tbl.stage(ind), 'late')
                    tidy_tbl.color(ind,:) = colors(5,:);
                end

                ind = ind+1;
            end
        end
    end
end

function [stat, p] = plot_width(avg)
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

    % Find the unique conditions in data
    genes = unique(avg.gene, 'stable');
%     reorder_ind = [2,4,3,1];
%     genes = genes(reorder_ind);

    stat = struct('p', cell(size(genes,1)+1,1),...
                  'tbl', [],...
                  'stats', [],...
                  'p_indiv', [],...
                  'means', [],... 
                  'names', []);

    p = cell(size(genes,1)+1,1);

    if size(genes,1) > 1
        f = figure('Position',[200,300,1000,300]);
        tiledlayout(f,1,size(genes,1));
    else
        figure
    end

    for i = 1:size(genes,1)
        % Find the unique conditions and initialize variables
        unordered_cond = unique(avg.condition(avg.gene == genes(i)), 'stable');
        cond_names = strings(size(unordered_cond));

        if strcmp(genes(i), 'zen')
            [stat(i), p{i}] = statistical_analysis(avg.gene_width(avg.gene == genes(i) & avg.stage == 'early'),...
                avg.condition(avg.gene == genes(i) & avg.stage == 'early'));
            [stat(end), p{end}] = statistical_analysis(avg.gene_width(avg.gene == genes(i) & avg.stage == 'late'),...
                avg.condition(avg.gene == genes(i) & avg.stage == 'late'));
            mean_x = [1:size(unordered_cond, 1);1:size(unordered_cond, 1)]';
            mean_y = zeros(size(unordered_cond, 1), 2);
            err_y = zeros(size(unordered_cond, 1), 2);
        else
            [stat(i), p{i}] = statistical_analysis(avg.gene_width(avg.gene == genes(i)), avg.condition(avg.gene == genes(i)));
            mean_x = (1:size(unordered_cond, 1))';
            mean_y = zeros(size(unordered_cond, 1), 1);
            err_y = zeros(size(unordered_cond, 1), 1);
        end
    
        % For each condition
        for j = 1:size(unordered_cond, 1)
            if strcmp(unordered_cond(j),'Control')
                cond_names(j) = unordered_cond(j);
            else
                cond_names(j) = append('\it',unordered_cond(j));
            end
            
            if strcmp(genes(i),'zen')
                mean_y(j,1) = mean(avg.gene_width((avg.gene == genes(i) & avg.condition == unordered_cond(j) & avg.stage == 'early')));
                err_y(j,1) = calc_error(avg.gene_width((avg.gene == genes(i) & avg.condition == unordered_cond(j) & avg.stage == 'early')), 'SD');
                mean_y(j,2) = mean(avg.gene_width((avg.gene == genes(i) & avg.condition == unordered_cond(j) & avg.stage == 'late')));
                err_y(j,2) = calc_error(avg.gene_width((avg.gene == genes(i) & avg.condition == unordered_cond(j) & avg.stage == 'late')), 'SD');
            else
                mean_y(j) = mean(avg.gene_width((avg.gene == genes(i) & avg.condition == unordered_cond(j))));
                err_y(j) = calc_error(avg.gene_width((avg.gene == genes(i) & avg.condition == unordered_cond(j))), 'SD');
            end
        end

        % Make a figure and keep the axis for plotting the raw data, the mean,
        % and the error bars on the same plot
        nexttile;
        hold on
        
        %colors(reorder_ind(i),:)
        % Plot raw data with jitter to offset the points with some transparency
        scatter(avg.x(avg.gene == genes(i)), avg.gene_width(avg.gene == genes(i)),...
            100, avg.color(avg.gene == genes(i),:), 'filled', 'jitter', 'on', 'jitteramount', 0.2,...
                'MarkerFaceAlpha', .5,'MarkerEdgeAlpha', .5);
        
        % Plot the mean with the error bars and set properties
        h = errorbar(mean_x, mean_y, err_y, '.', 'Color', 'k');
        set(h, 'linewidth', 2, 'markersize', 25);
        hold off

        set(gca, 'FontName', 'Arial')
        
        % Set properties of axis
        set(gca, 'xlim', [0, size(unordered_cond,1)+1],...
                 'xtick', 1:size(unordered_cond,1), ...
                 'xticklabels', cond_names,...
                 'XTickLabelRotation', 45,...
                 'ylim', [-0.025, 0.825],...%[0.11, 0.23] or [-0.025, 0.825]
                 'fontsize', 20);

        title(append('\it',genes(i)), 'FontName', 'Arial');
    end
end

function err_d = calc_error(d, select_error)
%CALC_ERROR Calculates the error for data d
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
%     STD_d = nanstd(d, [], 1);
    STD_d = std(d, 0, 1, 'omitnan');

    % standard error of the mean for data in d
    SEM_d = STD_d ./ sqrt(sum(~isnan(d), 1));

    % confidence interval for data in d
    ts_d = tinv(0.975, sum(~isnan(d), 1) - 1);
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
    
    % To calculate individual t tests between groups of two samples,
    % organize samples so samples to be tested are concatanated together
    % (for example column 1 and 2 will be tested, 3 and 4, etc)
    
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