function [x, y, g, c, groups, stat, p] = make_boxplot(field, gene, group_names, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~iscell(field)
    field = {field};
end

if ~iscell(gene)
    gene = {gene};
end

colors = [[0.12156862745098039, 0.4666666666666667, 0.7058823529411765];
          [1.0, 0.4980392156862745, 0.054901960784313725];
          [0.17254901960784313, 0.6274509803921569, 0.17254901960784313];
          [0.8392156862745098, 0.15294117647058825, 0.1568627450980392];
          [0.5803921568627451, 0.403921568627451, 0.7411764705882353];
          [0.5490196078431373, 0.33725490196078434, 0.29411764705882354];
          [0.8901960784313725, 0.4666666666666667, 0.7607843137254902];
          [0.4980392156862745, 0.4980392156862745, 0.4980392156862745];
          [0.7372549019607844, 0.7411764705882353, 0.13333333333333333];
          [0.09019607843137255, 0.7450980392156863, 0.8117647058823529]];

nuc_cols = {'A', 'B', 'M', 'mu', 'sig', 'dA', 'dB', 'dM', 'dmu', 'dsig'};
gene_cols = {'sV', 'sD', 'w', 'dsV', 'dsD', 'dw'};

data = cat(1, varargin{:});
nucproteins = struct('genotype', [], 'nucprotein_names', [], 'A', [],...
                     'B', [], 'M', [], 'mu', [], 'sig', [], 'dA', [],...
                     'dB', [], 'dM', [], 'dmu', [],...
                     'dsig', cell(16 .* size(data,1),1));
genes = struct('genotype', [], 'gene_names', [], 'sV', [], 'sD', [],...
               'w', [], 'dsV', [], 'dsD', [],...
               'dw', cell(16 .* size(data,1),1));
np = 1;
gn = 1;

for k = 1:size(data,1)
    for j = 1:size(data(k).nucprotein_names, 2)
        nucproteins(np).genotype = data(k).genotype;
        nucproteins(np).nucprotein_names = data(k).nucprotein_names{j};
        
        for i = 1:size(nuc_cols,2)
            nucproteins(np).(nuc_cols{i}) = data(k).(nuc_cols{i})(j);
        end
        
        np = np + 1;
    end
    
    for j = 1:size(data(k).gene_names, 2)
        genes(gn).genotype = data(k).genotype;
        genes(gn).gene_names = data(k).gene_names{j};
        
        for i = 1:size(gene_cols,2)
            genes(gn).(gene_cols{i}) = data(k).(gene_cols{i})(j);
        end
        
        gn = gn + 1;
    end
end

if np <= size(nucproteins, 1) && isempty(nucproteins(np).genotype)
    nucproteins(np:end) = [];
end

if gn <= size(genes, 1) && isempty(genes(gn).genotype)
    genes(gn:end) = [];
end

stat = cell(size(field,2));
p = cell(size(field,2));

for j = 1:size(field, 2)
    if any(strcmp(field{j}, nuc_cols))
        names = {nucproteins.nucprotein_names};
        val = [nucproteins.(field{j})];
        genotype = {nucproteins.genotype};
    elseif any(strcmp(field{j}, gene_cols))
        names = {genes.gene_names};
        val = [genes.(field{j})];
        genotype = {genes.genotype};
    else
        fprintf('%s not found.\n', field{j}')
        continue
    end

    stat{j} = cell(size(gene,2));
    p{j} = cell(size(gene,2));
    
    for i = 1:size(gene, 2)
        n = strcmp(names, gene{i});
        
        if any(n(:))
            y = val(n);
            g = genotype(n);
            
            groups = unique(g, 'stable');
            x = zeros(1, size(g,2));
            c = zeros(size(g,2), 3);
            mean_x = 1:size(groups, 2);
            mean_y = zeros(size(groups, 2),1);
            err_y = zeros(size(groups, 2),1);
            for k = 1:size(groups, 2)
                ind = strcmp(g, groups{k});
                x = k .* ind + x;
                c(ind,:) = bsxfun(@plus, c(ind,:), colors(mod(k-1, size(colors,1))+1,:));
                mean_y(k) = mean(y(ind));
                err_y(k) = calc_error(y(ind)','SD');
            end
            
            figure;
%             h = boxplot(y, g, 'width', .75, 'colors', [0 0 0],...
%                 'symbol', '', 'outliersize', 25);
            hold on
            scatter(x, y, 100, c, 'filled',...
                    'jitter', 'on', 'jitteramount', 0.3,...
                    'MarkerFaceAlpha', .5,'MarkerEdgeAlpha', .5);

             % Plot the mean with the error bars and set properties
            h = errorbar(mean_x, mean_y, err_y, '.', 'Color', 'k');
            set(h, 'linewidth', 2, 'markersize', 25);
            hold off
            title(sprintf('%s for %s', field{j}, gene{i}));
            set(h, 'linewidth', 2);

            if min(y) >= 0
                set(gca, 'fontsize', 20, 'ylim', [0, 1.05 .* max(y)]);
            else
                set(gca, 'fontsize', 20,...
                         'ylim', [1.05 .* min(y), 1.05 .* max(y)]);
            end

            set(gca, 'FontName', 'Arial')

            if isempty(group_names)
                group_names = groups;
            end
        
            % Set properties of axis
            set(gca, 'xlim', [0, size(groups,2)+1],...
                     'xtick', 1:size(groups,2), ...
                     'xticklabels', group_names,...
                     'XTickLabelRotation', 45,...
                     'fontsize', 20);
        else
            fprintf('%s not found for %s.\n', gene{i}, field{j})
        end

        [stat{j}{i}, p{j}{i}] = statistical_analysis(y, g);
    end
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