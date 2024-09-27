function [avg, tbl, stat, p, data] = stats_ms2_series(data, varargin)
%STATS_MS2_SERIES Plots the data from quantify_ms2
% 
%   Input
%       data: the output from quantify_ms2
%   Optional Inputs
%       varargin
%           1): condition, either a cell array of names or an array of the
%               indices of conditions in handle_inputs
%           2): two entry array corresponding to the x-axis limits
%           3): two entry array corresponding to the y-axis limits
%           4): enter a 1 or 2 to set the color scheme. 1 is chosen by
%               default
%           5): enter 'norm' to normalize data by the first time point, or
%               'unnorm' to plot unnormalized data. The default is to plot
%               the normalized data.
%           6): enter 'indiv' to specify plotting individual lines. The
%               default is to plot the average
%           7): enter the index that corresponds to the threshold entered
%               into quantify_ms2 (typically 1, 2, or 3). The default is 2.
%           8): enter 'SEM', 'CI', or 'SD' to plot error bars calculated as
%               standard error of the mean, confidence intervals, or
%               standard deviation. The default is SEM.
%           9): enter 'area plot' to plot area of expression
%           
% 
%   Output
%       None
% 
%   Overview
%       Plots the data from quantify_ms2, including number of spots,
%       average number of pixels per spot, and average mean intensity per
%       spot. The first input is the output from quantify_ms2. The other
%       inputs are optional. If no condition is given, the function will
%       try to plot all the conditions in handle_inputs. To add additional
%       conditions, add conditions and condition names, with colors, to
%       handle_inputs. data.blue_light is manually entered into the
%       structure as [z1,z2;t1,t2] where (z1,t1) is the z frame and time
%       frame where light is turned on, and (z2,t2) is the z frame and time
%       frame with the last detectable blue light. If no light is used,
%       blue_light should be NaN. 
% 
%   Note
%       The eighth input to varargin is not well tested for the rectangle,
%       and may result in errors. It is best to use the bar.
    
    % Check the entered condition against list of all conditions
    [condition, m, T] = handle_inputs(varargin{:});
    
    % Compute the averages
    [avg, tbl, data] = calc_means(condition, data, m, T, varargin{:});
    
%     % Perform ANOVA for determining statistical significance
%     if 0%size(avg, 1) > 1
% %         [stat, p] = statistical_analysis(avg);
%     else
%         stat = [];
%         p = [];
%     end

    plot_ms2(avg, varargin{:});

    % If plot boxplots is true
    if nargin >= 10 && strcmp(varargin{9}, 'area plot')
%         plot_boxplots(avg);
        [stat, p] = plot_domain(avg, varargin{:});
    else
        stat = [];
        p = [];
    end
end

function [condition, m, T] = handle_inputs(varargin)
%HANDLE_INPUTS Handles inputs for the conditions to be plotted
% 
%   Input
%       varargin: same varargin as entered to main function
%           1): entered conditions to plot data for
% 
%   Output
%       condition: updated condition with color and name specified
%       m: 1 or 2 to index normalized or unnormalized data respectively
%       T: index to select data for threshold (usually 1, 2 or 3)
% 
%   Overview
%       This function matches the entered condition with a condition name
%       and color. This allows manipulating the name that appears in the
%       legend during plotting, and the colors associated with a certain
%       condition. To add a condition, simply add a string to the array,
%       followed by a string of the name as it should be included in the
%       legend, followed by two color indicies that match the colors in
%       initiate_plot, separated by commas. Each condition should be
%       separated by semicolons. The second color is used by setting the
%       color scheme to 2 in the main function. This function also
%       determines if normalized or unnormalized data should be used, and
%       which threshold to plot based on user input to main function.
    
    % An array containing the condition, condition name, and two colors
    condition = {'wt early', 'sna.wt', 3, 3, 1;...
                 'wt late', 'sna.wt', 3, 3, 2;...
                 'ddis early', 'sna.\Deltadis', 2, 2, 3;
                 'ddis late', 'sna.\Deltadis', 2, 2, 4;...
                 'dprox early', 'sna.\Deltaprox', 4, 4, 5;...
                 'dprox late', 'sna.\Deltaprox', 4, 4, 6;...
                 'wt', 'wt', 3, 4, 1;...
                 'wt dprox', '\Deltaprox', 3, 5, 2;...
                 'wt ddis', '\Deltadis', 2, 2, 3;...
                 'dl1/cyo', 'dl^1/CyO', 6, 6, 4;...
                 'dl4/cyo', 'dl^4/CyO', 2, 4, 5;...
                 'dl4/cyo dprox', 'dl^4 \Deltaprox' 2, 5, 4;...
                 '30 sec', 'sna.wt 30 sec', 1, 1, 1;...
                 '1 min', 'dl-BLID/+ sna.wt 1 min', 4, 4, 2;...
                 'dark 3 min', 'dark 3 min', 5, 5, 3;...
                 '5 min late', '5 min late', 4, 4, 4;...
                 'dark', 'dark', 2, 4, 1;...%'dark', 'dl-BLID/+ sna.wt dark', 2, 4;...
                 '3 min', '3 min', 1, 4, 2;...%'3 min', 'dl-BLID/+ sna.wt 3 min', 1, 4;...
                 '5 min', '5 min', 3, 4, 3;...%'5 min', 'dl-BLID/+ sna.wt 5 min', 3, 4;...
                 'dprox dark', 'dark', 2, 5, 4;...%'dprox dark', 'dl-BLID/+ sna.\Deltaprox dark', 2, 5;...
                 'dprox 3 min', '3 min', 1, 5, 5;...%'dprox 3 min', 'dl-BLID/+ sna.\Deltaprox 3 min', 1, 5;...
                 'dprox 5 min', '5 min', 3, 5, 6;...%'dprox 5 min', 'dl-BLID/+ sna.\Deltaprox 5 min', 3, 5;...
                 'wt dark', 'wt dark', 2, 4, 7;...%'wt dark', 'wt sna.wt dark', 2, 4;...
                 'wt 5 min', 'wt 5 min', 3, 4, 8;...'wt 5 min', 'wt sna.wt 5 min', 3, 4;...
                 'wt dprox dark', '\Deltaprox dark', 4, 5, 3;...
                 'LEXY', 'LEXY', 1, 1, 1;...
                 'test', '\itdl-AAAA-Venus', 1, 1, 1;...
                 'test1', '\itdl-AAAA-Venus', 3, 1, 1;...
                 'control', '\itdl-Venus', 2, 2, 2;...
                 'control1', '\itdl-Venus', 4, 2, 2;...
                 'dl-LEXY dark', 'Dark', 8, 8, 1;...%2, 4, 1;
                 'dl-LEXY 10 min', 'Light 10 min', 1, 4, 2;...%4, 4, 2;
                 'dl-LEXY 20 min', 'Light 20 min', 4, 4, 3;...%3, 4, 3;
                 'dl-LEXY NC11 12', 'Light nc11-12', 9, 4, 3;%1, 4, 3;
                 'dl-LEXY NC11 13', 'Light nc11-13', 6, 4, 3;...%5, 4, 3;
                 'dl-LEXY NC13', 'Light nc13', 1, 1, 2;...%4, 4, 2;
                 'dl-LEXY NC13 weird', 'Light nc13', 4, 4, 2;...
                 'dl-LEXY NC12 14', 'Light nc12&14', 6, 6, 3;...%3, 4, 3;
                 'dl-LEXY NC12 14 weird', 'Light nc12&14', 3, 4, 3;...
                 'dl-LEXY NC13 14', 'Light nc13-14', 1, 4, 2;...%4, 4, 2;...
                 'dl-BLID dark', 'Dark', 8, 4, 1;...%2, 4, 1;
                 'dl-BLID 10 min', 'Light 10 min', 1, 4, 2;...%4, 4, 2;
                 'dl-BLID 20 min', 'Light 20 min', 4, 4, 3;...%3, 4, 3;
                 'dl-BLID NC13 14', 'Light nc13-14', 1, 4, 2;...
                 'dl-LEXY het dark', 'Dark', 2, 4, 1;...
                 'dl-LEXY het two light', 'Two Light Exposures', 4, 4, 2;...
                 'dl-LEXY het continuous', 'Continuous Light', 3, 4, 3;...
                 'dl-mCh-LEXY dark', 'Dark', 8, 4, 1;...%2, 4, 1;
                 'dl-mCh-LEXY NC13', 'Light nc13', 1, 4, 2;...%4, 4, 2;
%                  'dl-Venus', 'dl-Venus', 2, 4, 1;...
%                  'dl-AAAA-Venus', 'dl-AAAA-Venus', 4, 4, 2;...
                 'dl-LEXY dprox dark', 'Dark', 2, 4, 1;...
                 'dl-LEXY dprox 10 min', 'Light 10 min', 4, 4, 2;...
                 'dl-LEXY dprox 20 min', 'Light 20 min', 3, 4, 3;...
                 'dl-Venus', 'dl-Venus', 1, 5, 1;...
                 'dl-Venus/df', 'dl-Venus/df', 2, 2, 1;...
                 'dl-AAAA', 'dl-AAAA', 3, 3, 2;...
                 'dl-AAAA/df', 'dl-AAAA/df', 4, 4, 1;...
                 'sna p1.3', '\itsna \Deltaprox 1.3', 1, 3, 2;...
                 'sna d2.0', '\itsna \Deltadist 2.0', 2, 3, 2;...
                 'sna d0.2', '\itsna \Deltadist 0.2', 3, 3, 2};
    
    % If a condition was entered
    if nargin >= 1 && ~isempty(varargin{1})
        % If the condition is an index of the condition array
        if isnumeric(varargin{1})
            % Save all the settings for the entered condition
            condition = condition(varargin{1},:);
        % If the conditon is an array of strings
        elseif iscell(varargin{1})
            % Initialize looped variable
            index = zeros(size(condition,1),size(varargin{1},2));
            
            % For each entered condition
            for k = 1:size(varargin{1},2)
                % Compare the entered condition to the condition array
                index(:,k) = strcmp(varargin{1}{k}, condition(:,1));
            end
            
            % Save all the settings for the entered condition
            condition = condition(logical(sum(index,2)),:);
        end
    end
    
    % If unnorm was entered
    if nargin >= 5 && strcmp(varargin{5}, 'unnorm')
        % Set the index to pick unnormalized data
        m = 0;
    else
        % Otherwise set the index to pick normalized data
        m = 1;
    end
    
    % If a threshold index was given as an input
    if nargin >= 7 && ~isempty(varargin{7})
        % Set threshold index to entered index
        T = varargin{7};
    else
        % Otherwise set the index to second threshold
        T = 2;
    end
end

function [avg, tbl_tidy, data] = calc_means(condition, data, m, T, varargin)
%CALC_MEANS Calculates the averages for data
% 
%   Input
%       condition: condition returned from handle_inputs
%       data: data structure the user entered into main function
%       m: 1 or 2 to index normalized or unnormalized data respectively
%       T: index to select data for threshold (usualy 1, 2 or 3)
%       varargin: same varargin as entered into main function
% 
%   Output
%       avg: structure with averaged data
%       subset: subset of data structure containing only the data
%           matching the conditions
% 
%   Overview
%       This function calculates the average data for time, number of
%       spots, area, and intensity. In addition error bars are calculated
%       as the standard error of the mean, confidence intervals, or
%       standard deviation. The individual data for time, replicate number,
%       condition, area, and intensity per spot are also saved for
%       generating tidy data using export_ms2_quantification

    
%     % Initialize looped variables
%     n_reps = zeros(size(condition, 1), 1);
%     
%     % For each condition
%     for i = 1:size(condition, 1)
%         % save number of replicates
%         n_reps(i,1) = sum(strcmp({data.condition}, condition{i,1}));
%     end
    
%     if all(n_reps(:) == n_reps(1))
%         n_reps = n_reps(1);
%     end
    
%     n_t = size(cat(1, data.ind), 2);

%     field = {'t'; 'n'; 'A'; 'avg_I'; 'max_I'; 'sum_I'; 'group';...
%         'replicate'};
%     
%     initiate_single = repmat({zeros(n_t, 2, n_reps)}, 6, 1);
%     avg_single = table(initiate_single{:},...
%                        'VariableNames', field(1:6));
% 
%     initiate_pooled = repmat({zeros(n_t, 2)}, 6, 1);
%     avg_pooled = table(initiate_pooled{:},...
%                        'VariableNames', field(1:6));
%     
%     initiate_tidy = repmat({cell(n_t,1)}, 7, 1);
%     avg_tidy = table(initiate_tidy{:},...
%                      'VariableNames', field([1, 3:8]));
% 
%     avg = table(cell(size(condition, 1), size(condition, 2)),...
%                 zeros(size(condition, 1), 1),...
%                 repmat({zeros(2, 1)}, size(condition, 1), 1),...
%                 repmat({avg_single}, size(condition, 1), 1),...
%                 repmat({avg_pooled}, size(condition, 1), 1),...
%                 repmat({avg_tidy}, size(condition, 1), 1),...
%                 'VariableNames', {'condition',...
%                                   'total_n_reps',...
%                                   'blue',...
%                                   'single',...
%                                   'pooled',...
%                                   'tidy'});
    
    % Create a structure to hold avg data
    avg = struct('condition', cell(1,size(condition, 1)),...
                 'indiv', [], 'intrp', [], 'avg_intrp', [], ...
                 'avg_blue', NaN, 'blue', NaN,...
                 'max_n', [], 'avg_n', [], 'min_n', [], 'sum_n', [],...
                 'pts', [], 'hull_pts', [], 'hull_area', []);
    
    % Initialize variables
    j = false(size(condition, 1), size(data, 2));
    init_zeros = zeros(size(cat(1, data.time), 1), 1);
    init_string = strings(size(cat(1, data.time), 1), 1);
    init_tbl_tidy = {init_zeros, init_zeros, init_zeros, init_string};

    % Make a table using the cell array, so that there are t rows
    % and 2 column arrays in each table variable
    tbl_tidy = table(init_tbl_tidy{:},...
                    'VariableNames', {'t', 'n', 'replicate', 'condition'});
    
    % Index into the tidy data set
    tidy_ind = 1;

    % For each entered condition
    for i = 1:size(condition, 1)
        % Find indices of data that match entered condition (logical
        % indexing)
        j(i,:) = strcmp({data.condition}, condition{i,1});
        q = find(j(i,:));
        
        % Take all replicates that match condition
        reps = data(j(i,:));
        
        % Make cell arrays that are the size of the number of replicates
        avg(i).indiv = cell(size(reps,2), 1);
        avg(i).intrp = cell(size(reps,2), 1);
        avg(i).pts = cell(size(reps,2), 3);
        avg(i).hull_pts = cell(size(reps,2), 3);
        avg(i).ind = cell(size(reps,2), 3);

        % Initialize a cell array with 3 arrays containing zeros that
        % has the length of t rows and 2 columns
        init_tbl_nc = repmat({zeros(size(reps, 2), 1)}, 3, 1);

        % Make a table using the cell array, so that there are t rows
        % and 2 column arrays in each table variable
        avg(i).max_n = table(init_tbl_nc{:},...
                        'VariableNames', {'nc12', 'nc13', 'nc14',});
        avg(i).avg_n = table(init_tbl_nc{:},...
                        'VariableNames', {'nc12', 'nc13', 'nc14',});
        avg(i).min_n = table(init_tbl_nc{:},...
                        'VariableNames', {'nc12', 'nc13', 'nc14',});
        avg(i).sum_n = table(init_tbl_nc{:},...
                        'VariableNames', {'nc12', 'nc13', 'nc14',});
        avg(i).hull_area = table(init_tbl_nc{:},...
                        'VariableNames', {'nc12', 'nc13', 'nc14',});

        % Save the condition that was averaged
        avg(i).condition = condition(i,:);

        % For each replicate
        for r = 1:size(reps, 2)
%             avg(i).hull_pts{r} = cell(size(reps(r).time, 1), 1);
%             avg(i).hull_area{r} = zeros(size(reps(r).time, 1), 1);

            % Initialize a cell array with 6 arrays containing zeros that
            % has the length of t rows and 2 columns
            init_tbl_indiv = repmat({zeros(size(reps(r).time, 1), 2)}, 6, 1);

            % Make a table using the cell array, so that there are t rows
            % and 2 column arrays in each table variable
            avg(i).indiv{r} = table(init_tbl_indiv{:},...
                            'VariableNames', {'t', 'n', 'A', 'avg_I',...
                                              'max_I', 'sum_I'});

            % For each time point
            for t = 1:size(reps(r).time, 1)
                % Save the time
                avg(i).indiv{r}.t(t) = reps(r).time(t) - reps(r).time(reps(r).t_align);
                tbl_tidy.t(tidy_ind) = avg(i).indiv{r}.t(t);

                % Save the number of spots
                if isfield(reps, 'n_spots')
                    avg(i).indiv{r}.n(t) = reps(r).n_spots(t);
                    temp_n = size(reps(r).A{t,T}, 1);
                else
                    if isfield(reps, 'A')
                        avg(i).indiv{r}.n(t) = size(reps(r).A{t,T}, 1);
                        temp_n = size(reps(r).A{t,T}, 1);
                    end
                end
                tbl_tidy.n(tidy_ind) = avg(i).indiv{r}.n(t);
                
                % Save the replicate number and condition
                tbl_tidy.replicate(tidy_ind) = r;
                tbl_tidy.condition(tidy_ind) = avg(i).condition{1};

                % If a spot was detected
                if ~isempty(reps(r).avg_I{t}) || (avg(i).indiv{r}.n(t) ~= 0 && temp_n ~= 0)
                    % Calculate the average area of all spots at a certain
                    % time point, convert from number of pixels to
                    % microns^2. Calculate means for average intensity, max
                    % intensity, and sum of intensity across spots for each
                    % time point
                    if isfield(reps, 'A')
                        avg(i).indiv{r}.A(t,1) = mean(reps(r).A{t,T} .*...
                                            prod(reps(r).pixel_length,2));
                        avg(i).indiv{r}.A(t,2) = calc_error(reps(r).A{t,T} .*...
                                                prod(reps(r).pixel_length,2), varargin{8}, 1);
                    end

                    if isfield(reps, 'avg_I')
                        avg(i).indiv{r}.avg_I(t,1) = mean(reps(r).avg_I{t,T});
                        avg(i).indiv{r}.avg_I(t,2) = calc_error(reps(r).avg_I{t,T}, varargin{8}, 1);
                    end

                    if isfield(reps, 'max_I')
                        avg(i).indiv{r}.max_I(t,1) = mean(reps(r).max_I{t,T});
                        avg(i).indiv{r}.max_I(t,2) = calc_error(reps(r).max_I{t,T}, varargin{8}, 1);
                    end

                    if isfield(reps, 'sum_I')
                        avg(i).indiv{r}.sum_I(t,1) = mean(reps(r).sum_I{t,T});
                        avg(i).indiv{r}.sum_I(t,2) = calc_error(reps(r).sum_I{t,T}, varargin{8}, 1);
                    end

                    % Calculate the error determined by the inputed value
                    % for area of spots, average intensity, max intensity
                    % and sum of intensities
                    
                    
                end
                
                % Increment the tidy data set index
                tidy_ind = tidy_ind + 1;
            end

            % If normalized
            if m == 1 %&& ~isnan(reps(r).t_norm)
%                 time_norm = 'max 12';
%                 time_norm = 'max 13';
                time_norm = 'input';
%                 max_field = 'n';
                max_field = 'avg_I';

                % If normalzing by value of inputed time
                if strcmp(time_norm, 'input')
                    % Save index to normalize data by
                    tn = reps(r).t_norm;
                % Elseif normalizing by maximum value
                elseif strcmp(time_norm, 'max 12')
                    % Save maximum value index to normalize data by
                    [~, tn] = max(avg(i).indiv{r}.(max_field)(reps(r).nuc_cycle(1,1):reps(r).nuc_cycle(1,2),1));
                    tn = tn + reps(r).nuc_cycle(1,1) - 1;
                % Elseif normalizing by maximum value
                elseif strcmp(time_norm, 'max 13')
                    % Save maximum value index to normalize data by
                    [~, tn] = max(avg(i).indiv{r}.(max_field)(reps(r).nuc_cycle(2,1):reps(r).nuc_cycle(2,2),1));
                    tn = tn + reps(r).nuc_cycle(2,1) - 1;
                end

                % Normalize data by value saved in index tn
                if isfield(reps, 'A')
                    avg(i).indiv{r}.n = avg(i).indiv{r}.n ./...
                                        avg(i).indiv{r}.n(tn,1);
                    avg(i).indiv{r}.A = avg(i).indiv{r}.A ./...
                                        avg(i).indiv{r}.A(tn,1);
                end

                if isfield(reps, 'avg_I')
                    avg(i).indiv{r}.avg_I = avg(i).indiv{r}.avg_I ./...
                                            avg(i).indiv{r}.avg_I(tn,1);
                end

                if isfield(reps, 'max_I')
                    avg(i).indiv{r}.max_I = avg(i).indiv{r}.max_I ./...
                                            avg(i).indiv{r}.max_I(tn,1);
                end

                if isfield(reps, 'max_I')
                    avg(i).indiv{r}.sum_I = avg(i).indiv{r}.sum_I ./...
                                            avg(i).indiv{r}.sum_I(tn,1);
                end
            end
            
            if nargin >= 13 && strcmp(varargin{9}, 'area plot')
                data(q(r)).rm_pts = data(q(r)).centers;
                % For 3 nuclear cycles
                for k = 1:3
                    % save range of time poiint indices using inputed nuclear
                    % cycle start and end indices
                    t_range = reps(r).nuc_cycle(k,1):reps(r).nuc_cycle(k,2);
    
                    % If nc index is not NaN
                    if ~isnan(reps(r).nuc_cycle(k,1)) && ~isnan(reps(r).nuc_cycle(k,2))
                        % Find max, mean, min, and area under the curve for the
                        % given time range defined by the nc index
                        avg(i).max_n{r,k} = max(avg(i).indiv{r}.n(t_range,1));
                        avg(i).avg_n{r,k} = mean(avg(i).indiv{r}.n(t_range,1));
                        avg(i).min_n{r,k} = min(avg(i).indiv{r}.n(t_range,1));
                        avg(i).sum_n{r,k} = trapz(avg(i).indiv{r}.t(t_range,1),avg(i).indiv{r}.n(t_range,1));

                        [avg(i).pts{r,k}, avg(i).ind{r,k}] = rmoutliers(cat(1, reps(r).centers{t_range}), 'ThresholdFactor', 2);
%                         avg(i).pts{r,k} = rmoutliers(cat(1, reps(r).centers{t_range}), 'ThresholdFactor', 3);
                        [avg(i).hull_pts{r,k},avg(i).hull_area{r,k}] = convhull(avg(i).pts{r,k});
                        avg(i).hull_area{r,k} = avg(i).hull_area{r,k} .* prod(reps(r).pixel_length);
                        
%                         data(q(r)).pts{k} = avg(i).pts{r,k};
%                         data(q(r)).hull_pts{k} = avg(i).hull_pts{r,k};
                        p = 1;

                        for qqq = t_range
                            if ~isempty(data(q(r)).rm_pts{qqq})
                               data(q(r)).temp{qqq} = ~avg(i).ind{r,k}(p:(p+size(data(q(r)).rm_pts{qqq}, 1) - 1));
                                data(q(r)).rm_pts{qqq} = data(q(r)).rm_pts{qqq}(~avg(i).ind{r,k}(p:(p+size(data(q(r)).centers{qqq}, 1) - 1)),:);
                                p = p + size(data(q(r)).centers{qqq}, 1);
                            end
                        end

                    % Else, if nc index is NaN
                    else
                        % Set value to NaN
                        avg(i).max_n{r,k} = NaN;
                        avg(i).avg_n{r,k} = NaN;
                        avg(i).min_n{r,k} = NaN;
                        avg(i).sum_n{r,k} = NaN;
                    end
                end

                data(q(r)).pts = avg(i).pts(r,:);
                data(q(r)).hull_pts = avg(i).hull_pts(r,:);
            end
        end

        % Find min and max time
        min_t = min(tbl_tidy.t(strcmp(tbl_tidy.condition, avg(i).condition{1})));
        max_t = max(tbl_tidy.t(strcmp(tbl_tidy.condition, avg(i).condition{1})));
        
        intrp_field = {'n','A','avg_I','max_I','sum_I'};
%             intrp_field = {'avg_I'};
        % Initilize variables
        temp_intrp = zeros(size((min_t:0.2:max_t)', 1), size(reps, 2), size(intrp_field,2));

        % Initialize table by replicating an array of zeros for each
        % variable
        init_tbl_intrp = repmat({zeros(size((min_t:0.2:max_t)', 1), 1)}, 6, 1);

        avg(i).avg_intrp = table(init_tbl_intrp{:},...
                            'VariableNames', {'t', 'n', 'A', 'avg_I',...
                                              'max_I', 'sum_I'});
        
        % For each replicate
        for r = 1:size(reps, 2)
            % Make a table using the cell array, so that there are t rows
            % and 2 column arrays in each table variable
            avg(i).intrp{r} = table(init_tbl_intrp{:},...
                            'VariableNames', {'t', 'n', 'A', 'avg_I',...
                                              'max_I', 'sum_I'});
            
            % Save the interperolated time
            avg(i).intrp{r}.t = (min_t:0.2:max_t)';

            for w = 1:size(intrp_field,2)
                % Interpolate the data for the given time
                avg(i).intrp{r}.(intrp_field{w}) = interp1(avg(i).indiv{r}.t(:,1),...
                    avg(i).indiv{r}.(intrp_field{w})(:,1), avg(i).intrp{r}.t, 'makima', NaN);
                % Save the interpolated values in an array for averaging
                temp_intrp(:,r,w) = avg(i).intrp{r}.(intrp_field{w});
            end
        end

        for w = 1:size(intrp_field,2)
            % Conctenate the time, mean, and error for the interpolated values
            avg(i).avg_intrp.t = (min_t:0.2:max_t)';
            avg(i).avg_intrp.(intrp_field{w}) = cat(2, mean(temp_intrp(:,:,w), 2, 'omitnan'), calc_error(temp_intrp(:,:,w), varargin{8}, 2));
        end

        
%         % For each time point
%         for x = 1:size(ind, 1)
%             % Initialize variables that are only saved for replicates
%             A = cell(size(ind, 2), 1);
%             avg_I = cell(size(ind, 2), 1);
%             max_I = cell(size(ind, 2), 1);
%             sum_I = cell(size(ind, 2), 1);
%             indiv_t = cell(size(ind, 2), 1);
%             replicate = cell(size(ind, 2), 1);
%             
%             % For each replicate
%             for y = 1:size(ind, 2)
%                 % Store data in the range set by ind. Subtract the time in
%                 % the first index in ind to set beginning of plot to zero
%                 avg.single{i}.t(x,1,y) = reps(y).time(ind(x,y)) -...
%                                                 reps(y).time(ind(1,y));
%                 
%                 % Store data on number of MS2 spots in the range set by
%                 % ind, for normalized or unnormalized data (m) and for the
%                 % threshold (T)
%                 avg.single{i}.n(x,1,y) = reps(y).n_spots(ind(x,y), T);
%                 
%                 % Store individual MS2 spot areas in the range set by ind,
%                 % for the threshold (T), converting area from number of
%                 % pixels to the pixel length usually in microns.
%                 A{y} = reps(y).A{ind(x,y), T} .*...
%                         prod(reps(y).pixel_length,2);
% 
%                 % Store individual MS2 spot mean intensity, max intensity,
%                 % and sum of intensities in the range set by ind, for the
%                 % threshold (T)
%                 avg_I(y) = reps(y).avg_I(ind(x,y), T);
%                 max_I(y) = reps(y).max_I(ind(x,y), T);
%                 sum_I(y) = reps(y).sum_I(ind(x,y), T);
%                 
%                 % If any spots were detected
%                 if avg.single{i}.n(x,1,y) ~= 0 
%                     % Find the mean for area, mean intesity, max intensity,
%                     % and sum of intensities
%                     avg.single{i}.A(x,1,y) = mean(A{y});
%                     avg.single{i}.avg_I(x,1,y) = mean(avg_I{y});
%                     avg.single{i}.max_I(x,1,y) = mean(max_I{y});
%                     avg.single{i}.sum_I(x,1,y) = mean(sum_I{y});
%                     
%                     % Compute error for means
%                     avg.single{i}.A(x,2,y) = calc_error(A{y}, varargin{9});
%                     avg.single{i}.avg_I(x,2,y) = calc_error(avg_I{y},...
%                                                             varargin{9});
%                     avg.single{i}.max_I(x,2,y) = calc_error(max_I{y},...
%                                                             varargin{9});
%                     avg.single{i}.sum_I(x,2,y) = calc_error(sum_I{y},...
%                                                             varargin{9});
%                 end
%                 
%                 % Assign the time for each image to each spot in that image
%                 indiv_t{y} = avg.single{i}.t(x,1,y) .* ones(size(A{y}));
%                 
%                 % Assign a replicate number to each spot
%                 replicate{y} = y .* ones(size(A{y}));
%             end
% 
%             if ~all(avg.single{i}.n(x,1,:) == 0)
%                 % Combine all individual spot times, areas, mean
%                 % intensities, max intensities, and sum intensities
%                 avg.tidy{i}.t{x,1} = cat(1, indiv_t{:});
%                 avg.tidy{i}.A{x,1} = cat(1, A{:});
%                 avg.tidy{i}.avg_I{x,1} = cat(1, avg_I{:});
%                 avg.tidy{i}.max_I{x,1} = cat(1, max_I{:});
%                 avg.tidy{i}.sum_I{x,1} = cat(1, sum_I{:});
%                 
%                 % Compute the average for number of pixels per focus,
%                 % mean intensity, max intensity, and sum of intensities
%                 avg.pooled{i}.A(x,1) = mean(avg.tidy{i}.A{x,1}, 1,...
%                                             'omitnan');
%                 avg.pooled{i}.avg_I(x,1) = mean(avg.tidy{i}.avg_I{x,1},...
%                                                 1, 'omitnan');
%                 avg.pooled{i}.max_I(x,1) = mean(avg.tidy{i}.max_I{x,1},...
%                                                 1, 'omitnan');
%                 avg.pooled{i}.sum_I(x,1) = mean(avg.tidy{i}.sum_I{x,1},...
%                                                 1, 'omitnan');
% 
%                 % Combine all individual spot replicate numbers
%                 avg.tidy{i}.replicate{x,1} = cat(1, replicate{:});
% 
%                 % Assign a condition to each spot
%                 avg.tidy{i}.group{x,1} = repmat(condition(i,1),...
%                                            size(avg.tidy{i}.A{x,1}));
% 
%                 % Compute the error bars
%                 avg.pooled{i}.A(x,2) = calc_error(avg.tidy{i}.A{x,1},...
%                                                   varargin{9});
%                 avg.pooled{i}.avg_I(x,2) = calc_error(...
%                                     avg.tidy{i}.avg_I{x,1}, varargin{9});
%                 avg.pooled{i}.max_I(x,2) = calc_error(...
%                                     avg.tidy{i}.max_I{x,1}, varargin{9});
%                 avg.pooled{i}.sum_I(x,2) = calc_error(...
%                                     avg.tidy{i}.sum_I{x,1}, varargin{9});
% 
%                 % If data should be normalized
%                 if m == 1
%                     % Normalize data by the mean of the first time point
%                     avg.tidy{i}.A{x,1} = avg.tidy{i}.A{x,1} ./...
%                                             avg.pooled{i}.A(1,1);
%                     avg.tidy{i}.avg_I{x,1} = avg.tidy{i}.avg_I{x,1} ./...
%                                                 avg.pooled{i}.avg_I(1,1);
%                     avg.tidy{i}.max_I{x,1} = avg.tidy{i}.max_I{x,1} ./...
%                                                 avg.pooled{i}.max_I(1,1);
%                     avg.tidy{i}.sum_I{x,1} = avg.tidy{i}.sum_I{x,1} ./...
%                                                 avg.pooled{i}.sum_I(1,1);
%                 end
%             end
%         end
%         
%         % Save the condition that was averaged
%         avg.condition(i,:) = condition(i,:);
%         
%         % Compute the average for time and number of spots
%         avg.pooled{i}.t(:,1) = mean(squeeze(avg.single{i}.t(:,1,:)), 2);
%         avg.pooled{i}.n(:,1) = mean(squeeze(avg.single{i}.n(:,1,:)), 2,...
%                                     'omitnan');
%         
%         % Compute the error of the mean for time and number of spots
%         avg.pooled{i}.t(:,2) = calc_error(squeeze(...
%                                 avg.single{i}.t(:,1,:))', varargin{9})';
%         avg.pooled{i}.n(:,2) = calc_error(squeeze(...
%                                 avg.single{i}.n(:,1,:))', varargin{9})';
% 
%         % If normalized
%         if m == 1
%             % Normalize data by the value for the first time point
%             avg.single{i}.n = avg.single{i}.n ./ avg.single{i}.n(1,1,:);
%             avg.single{i}.A = avg.single{i}.A ./ avg.single{i}.A(1,1,:);
%             avg.single{i}.avg_I = avg.single{i}.avg_I ./...
%                                     avg.single{i}.avg_I(1,1,:);
%             avg.single{i}.max_I = avg.single{i}.max_I ./...
%                                     avg.single{i}.max_I(1,1,:);
%             avg.single{i}.sum_I = avg.single{i}.sum_I ./...
%                                     avg.single{i}.sum_I(1,1,:);
% 
%             avg.pooled{i}.n = avg.pooled{i}.n ./ avg.pooled{i}.n(1,1);
%             avg.pooled{i}.A = avg.pooled{i}.A ./ avg.pooled{i}.A(1,1);
%             avg.pooled{i}.avg_I = avg.pooled{i}.avg_I ./...
%                                     avg.pooled{i}.avg_I(1,1);
%             avg.pooled{i}.max_I = avg.pooled{i}.max_I ./...
%                                     avg.pooled{i}.max_I(1,1);
%             avg.pooled{i}.sum_I = avg.pooled{i}.sum_I ./...
%                                     avg.pooled{i}.sum_I(1,1);
%         end
%         
%         % Save the number of lines that were averaged (usually 3)
%         avg.total_n_reps(i) = size(avg.single{i}.t, 3);
        
%         % For plotting individual lines, save the condition information
%         % by overwriting condition in the data structure
%         [data(j(i,:)).condition] = deal(condition(i,:));
        
        % If no blue light illumination was used
        if all(isnan([data(j(i,:)).blue_light]))
            % Specify no blue light using NaN
            avg(i).avg_blue = [NaN; NaN];
        else
            % Otherwise average the blue light start and stop times for all
            % data points that have the given condition
            
            % Initilize variables
            blue_t = nan(2, size(reps, 2), size(reps(1).blue_light, 3));
            
            % For each replicate
            for r = 1:size(reps, 2)
                % Determine time of illumination
                % line 3) index using sub2ind, with the size of raw_time
                % line 4) first dimension is channel, and is always 1
                % line 5) second dimension is z-slice, the first row of
                %         blue_light
                % line 6) third dimension is time point, the 2nd row of
                %         blue_light
                % line 7) subtract 1 from the linear index, illumination
                %         begins at the end of the last time point
                % line 8) subtract the average time for the first index to
                %         make sure time is normalized properly
                for b = 1:size(reps(r).blue_light, 3)
                    blue_t(:,r,b) =  reps(r).raw_time(...
                                         (sub2ind(size(reps(r).raw_time),...
                                                  [1, 1],...
                                                  reps(r).blue_light(1,:,b),...
                                                  reps(r).blue_light(2,:,b))...
                                         - [1,0])'...
                                        ) - reps(r).time(reps(r).t_align);
                end
            end
            
            % Take the average blue light window
            avg(i).avg_blue = squeeze(mean(blue_t,2));
            avg(i).blue = blue_t;
        end
    end
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

% function [stat, p] = statistical_analysis(avg)
% %STATISTICAL_ANALYSIS Perform ANOVA to compare the means between conditions
% % 
% %   Input
% %       avg: the structure returned from calc_means
% % 
% %   Output
% %       stat: structure containing outputs from anova1 and multcompare
% %       pn: table of p-values for pairwise comparisons of number of spots
% %           at each timepoint
% %       pA: table of p-values for pairwise comparisons of areas at each
% %           timepoint
% %       pI: table of p-values for pairwise comparisons of intensities at
% %           each timepoint
% % 
% %   Overview
% %       This function performs statistical analysis on the data, number of
% %       spots, area, and intensity. Specifically, it performs one way ANOVA
% %       using anova1 and multiple comparisons using Tukey's HSD using
% %       multcompare. It returns the outputs from anova1 and multcompare in
% %       the structure stat and tables of p-values pA and pI for the
% %       pairwise comparison between conditions.
%     
%     f = {'n', 'A', 'avg_I', 'max_I', 'sum_I'};
% 
%     initiate_table = table('Size',[size(avg.tidy{1},1), sum(1:(size(avg, 1)-1))],...
%                 'VariableTypes', repmat({'double'}, 1, sum(1:(size(avg, 1)-1))));
%     p = struct('n', initiate_table, 'A', initiate_table, 'avg_I', initiate_table, 'max_I', initiate_table, 'sum_I', initiate_table);
% 
%     try
%         % Try to concatenate average data
%         condition = avg.condition;
%         singles = cat(1, avg.single{:});
%         n = permute(reshape(squeeze(singles.n(:,1,:)),...
%                             size(avg.tidy{1},1),...
%                             size(avg, 1),...
%                             avg.total_n_reps(1)), [3,2,1]);
% %         p.n = cell(size(n, 2), size(n, 2), size(n, 3));
% 
%         tidy_data = cat(1,avg.tidy{:});
% 
%         data = cell(size(avg.tidy{1},1), size(avg, 1), size(f, 2));
%         group = reshape(tidy_data.group, size(avg.tidy{1},1), size(avg, 1));
% 
%         for i = 2:size(f,2)
%             data(:,:,i) = reshape(tidy_data.(f{i}), size(avg.tidy{1},1), size(avg, 1));
%             
% %             p.(f{i}) = cell(size(data, 2), size(data, 2), size(data, 1));
% %             new_p.(f{i}) = table('Size',[size(data, 1), sum(1:(size(data, 2)-1))],...
% %                 'VariableTypes', repmat({'double'}, 1, sum(1:(size(data, 2)-1))));
% %             A = reshape(tidy_data.A, size(avg.tidy{1},1), size(avg, 1));
% %             I = reshape(tidy_data.avg_I, size(avg.tidy{1},1), size(avg, 1));
% %             max_I = reshape(tidy_data.max_I, size(avg.tidy{1},1), size(avg, 1));
% %             sum_I = reshape(tidy_data.sum_I, size(avg.tidy{1},1), size(avg, 1));
% %             group = reshape(tidy_data.group, size(avg.tidy{1},1), size(avg, 1));
%         end
%     catch
%         % If concatenation fails, the number of timepoints aren't equal
%         error('Number of time points are not equal');
%     end
% 
%     % Initialize a structure for storing the results of the statistical
%     % analysis
%     initiate_stat = repmat({cell(size(avg.tidy{1},1), 1)}, 6, 1);
%     stat_table = table(initiate_stat{:},...
%                        'VariableNames', {'p',...
%                                          'tbl',...
%                                          'stats',...
%                                          'p_indiv',...
%                                          'means',...
%                                          'names'});
%     struct_inputs = cat(1, f, ...
%                            repmat({stat_table}, 1, size(f, 2)));
%     stat = struct(struct_inputs{:});
% 
%     
% %     stat = struct('p_n', cell(size(A, 1), 1),...
% %                   'tbl_n', [],...
% %                   'stats_n', [],...
% %                   'p_indiv_n', [],...
% %                   'means_n', [],... 
% %                   'names_n', [],...
% %                   'p_A', [],...
% %                   'tbl_A', [],...
% %                   'stats_A', [],...
% %                   'p_indiv_A', [],...
% %                   'means_A', [],... 
% %                   'names_A', [],...
% %                   'p_I', [],...
% %                   'tbl_I', [],...
% %                   'stats_I', [],...
% %                   'p_indiv_I', [],...
% %                   'means_I', [],... 
% %                   'names_I', []);
% 
% %     p = struct('n', [], 'A', [], 'avg_I', [], 'max_I', [], 'sum_I', []);
% %     
% %     % Initialize variables for making comparison tables
% %     
% %     p.n = cell(size(n, 2), size(n, 2), size(n, 3));
% %     p.A = cell(size(A, 2), size(A, 2), size(A, 1));
% %     p.avg_I = cell(size(I, 2), size(I, 2), size(I, 1));
% %     p.max_A = cell(size(A, 2), size(A, 2), size(A, 1));
% %     p.sum_I = cell(size(I, 2), size(I, 2), size(I, 1));
%     
%     % For each time point
%     for i = 1:size(avg.tidy{1}, 1)
%         % Perform ANOVA on the number of spots grouped by condition
%         [stat.n.p{i}, stat.n.tbl{i}, stat.n.stats{i}, stat.n.p_indiv{i},...
%             stat.n.means{i}, stat.n.names{i},...
%             ] = stat_test(n(:,:,i,1), condition(:,1)');
%         
%         for j = 2:size(f, 2)
%              % Perform ANOVA on the area data grouped by condition
%             [stat.(f{j}).p{i}, stat.(f{j}).tbl{i}, stat.(f{j}).stats{i}, stat.(f{j}).p_indiv{i},...
%                 stat.(f{j}).means{i}, stat.(f{j}).names{i},...
%                 ] = stat_test(cat(1, data{i,:,j}), cat(1, group{i,:}));
%             
%             table_var_names = cell(sum(1:(size(data, 2)-1)),1);
%             % For each comparison
%             for k = 1:size(stat.(f{j}).p_indiv{i}, 1)
%                 % save the p-value for area in the p-value table
% %                 p.(f{j}){stat.(f{j}).p_indiv{i}(k,1) + 1,...
% %                     stat.(f{j}).p_indiv{i}(k,2), i} = stat.(f{j}).p_indiv{i}(k,6);
% 
%                 p.(f{j}){i,k} = stat.(f{j}).p_indiv{i}(k,6);
%                 table_var_names(k,1) = join({stat.(f{j}).names{i}{stat.(f{j}).p_indiv{i}(k,1)},...
%                     stat.(f{j}).names{i}{stat.(f{j}).p_indiv{i}(k,2)}}, ' compared to ');
%                 
%     %             % save the p-value for area in the p-value table
%     %             p.A{stat(i).p_indiv_A(k,1) + 1,...
%     %                 stat(i).p_indiv_A(k,2), i} = stat(i).p_indiv_A(k,6);
%     %             
%     %             % save the p-value for intensity in the p-value table
%     %             p.I{stat(i).p_indiv_I(k,1) + 1,...
%     %                 stat(i).p_indiv_I(k,2), i} = stat(i).p_indiv_I(k,6);
%             end
% 
%             p.(f{j}).Properties.VariableNames = table_var_names;
%         end
%         
% %         % Perform ANOVA on the area data grouped by condition
% %         [stat(i).p_A, stat(i).tbl_A, stat(i).stats_A, stat(i).p_indiv_A,...
% %             stat(i).means_A, stat(i).names_A,...
% %             p.A] = stat_test(cat(1, A{i,:}), cat(1, group{i,:}), p.A, i);
% %         
% %         % Perform ANOVA on the intensity data grouped by condition
% %         [stat(i).p_I, stat(i).tbl_I, stat(i).stats_I, stat(i).p_indiv_I,...
% %             stat(i).means_I, stat(i).names_I,...
% %             p.I] = stat_test(cat(1, I{i,:}), cat(1, group{i,:}), p.I, i);
%     end
% end
% 
% function [p, tbl, stats, p_indiv, means, names] = stat_test(data,...
%     group)
% %STAT_TEST Perform ANOVA to compare the means between conditions
% % 
% %   Input
% %       data: data that anova will be performed on
% %       group: identifier for data to correctly group it
% %       p_tbl: the p-value array for storing p-values of mutiple
% %              comparisons
% %       i: the index specifying which time point
% % 
% %   Output
% %       p: p-value from the anova
% %       tbl: a table returned from anova
% %       stats: statistics for mutiple comparison tests
% %       p_indiv: pairwise p-values from mutiple comparisons
% %       means: estimated means
% %       names: names of groups
% %       p_tbl: a table of p-values from mutiple comparisons
% % 
% %   Overview
% %       This function performs statistical analysis on the data. 
% %       Specifically, it performs one way ANOVA using anova1 and multiple
% %       comparisons using Tukey's HSD using multcompare. It returns the
% %       outputs from anova1 and multcompare in the structure stat and a
% %       table of p-values from the pairwise comparison between conditions.
% 
%     % Perform ANOVA on the data grouped by condition in group
%     [p, tbl, stats] = anova1(data, group, 'off');
%     
%     % Perform pairwise comparisons of data between conditions
%     % using Tukey's HSD
%     [p_indiv, means, ~, names] = multcompare(stats, 'display', 'off');
%     
% %     % Save time point i in table of mutiple comparisons
% %     p_tbl{1, 1, i} = sprintf('t%d', i);
% %     
% %     % Make row names of conditions for comparison
% %     p_tbl(2:end, 1, i) = names(1:(end-1));
% %     
% %     % Make column names of conditions for comparison
% %     p_tbl(1, 2:end, i) = names(2:end);
% end

function plot_ms2(avg, varargin)
%PLOT_AVG_MS2 Plots the average of the data.
% 
%   Input
%       avg: average structure computed in calc_means function
%       varargin: same varargin as entered into main function
% 
%   Output
%       None
% 
%   Overview
%       Plots the average of the data specified by the conditions.
%       Also plots a blue bar if blue light was used, based on the data
%       provided by data.blue_light. See main fuction description for how
%       blue_light data should be entered into the data structure. Instead
%       of a blue bar, a blue reactangle can also be plotted, however the
%       bar is preferred, and the rectangle may not work unless the data is
%       organized a certain way, due to lack of extensive testing.
    
    % Initiate the plot by defining array of colors, names of plots, names
    % of fields to plot, and determing where the blue bars should be
    % plotted
    [colors, x_label, y_label, y_blue, field] = initiate_plot(...
                                                4, varargin{:}); %size(avg, 2)
    
    % For each field
    for j = 1:size(field, 2)
        % Initilize variables in loop and make figure
        legend_names = cell(1, size(avg, 2));
        h_avg = zeros(size(avg, 2),1);
        h_indiv_legend = zeros(size(avg, 2),1);
        
        % If plotting average
        if isempty(varargin{6}) || any(strcmp(varargin{6}, 'avg')) || any(strcmp(varargin{6}, 'indiv'))
            % Make figure and use hold on to prevent overwriting plots
            figure;
            ax_avg = axes;
            hold(ax_avg, 'on');
        end

        % For each condition
        for i = 1:size(avg, 2)
            % If color scheme entered is 2
            if (nargin >= 5) && ~isempty(varargin{4}) && (varargin{4} >= 2)
                % Set color to color scheme to inputted number
                c = avg(i).condition{2+varargin{4}};
            else
                % Otherwise set color to color scheme 1
                c = avg(i).condition{3};
            end
            
            % Save the legend name based on the condition
            legend_names{i} = avg(i).condition{2};

%             if m == 1
%                 norm_test = mean(data(i).(field{j,1})(data(i).ind(1)), T);
%             else
%                 norm_test = 1;
%             end
            
            % If plotting average
            if isempty(varargin{6}) || any(strcmp(varargin{6}, 'avg'))
                % Save time, values/means, and error
                x = avg(i).avg_intrp.t;
                y = avg(i).avg_intrp.(field{j})(:,1);
                dy = avg(i).avg_intrp.(field{j})(:,2);
    
                % Create a filled object to display error as a shaded region
                p = fill(ax_avg, [x;flipud(x)], [y-dy;flipud(y+dy)], colors(c,:), 'linestyle', 'none');
    
                % Set transparency to 50%
                alpha(p, 0.15);
    
                % Plot mean values
                h_avg(i) = line(ax_avg, x, y, 'color', colors(c,:));
    
                % Set line width and marker size
                set(h_avg(i), 'linewidth', 4, 'markersize', 15);
            end

            % If plotting individual lines
            if ~isempty(varargin{6}) && (any(strcmp(varargin{6}, 'indiv')) || any(strcmp(varargin{6}, 'indiv_separate')))
                if ~isempty(varargin{6}) && any(strcmp(varargin{6}, 'indiv_separate'))
                    % Make figure and use hold on to prevent overwriting plots
                    figure;
                    ax_indiv = axes;
                    hold(ax_indiv, 'on');

                    % To generate lines of the same color but are brighter or
                    % darker, make beta, a factor for scaling
                    % brightness/darkness
                    beta = (-0.1 * (size(avg(i).indiv, 1) - 1)):0.2:(0.1 * (size(avg(i).indiv, 1) - 1));
                else
                    ax_indiv = ax_avg;
                    beta = [];
                end

                % Initialize variables
                h_indiv = zeros(size(avg(i).indiv, 1), 1);

                % For each replicate
                for k = 1:size(avg(i).indiv, 1)
                    % If beta is greater than zero
                    if ~isempty(beta) && beta(k) > 0
                        % Make colors lighter
                        new_color = colors(c,:).^(1-beta(k));
                    % If beta is less than or equal to 0
                    elseif ~isempty (beta) && beta(k) <= 0
                        % Make colors darker
                        new_color = colors(c,:).^(1/(1+beta(k)));
                    else
                        new_color = colors(c,:);
                    end
    
                    % Plot the average data with error bars, time versus data
                    % specified by field
                    
    %                 h3(k) = errorbar(avg(i).indiv{k}.t(:,1),...
    %                                  avg(i).indiv{k}.(field{j})(:,1),...
    %                                  avg(i).indiv{k}.(field{j})(:,2),...
    %                                  avg(i).indiv{k}.(field{j})(:,2),...
    %                                  avg(i).indiv{k}.t(:,2),...
    %                                  avg(i).indiv{k}.t(:,2),...
    %                                  '.-', 'Color', new_color);
                    
                    % Save time, values/means, and error
                    x = avg(i).indiv{k}.t(:,1);
                    y = avg(i).indiv{k}.(field{j})(:,1);
                    dy = avg(i).indiv{k}.(field{j})(:,2);

                    % Create a filled object to display error as a shaded region
                    p = fill(ax_indiv, [x;flipud(x)], [y-dy;flipud(y+dy)], new_color, 'linestyle', 'none');

                    % Set transparency to 50%
                    alpha(p, 0.15);

                    % Plot mean values
                    h_indiv(k) = line(ax_indiv, x, y, 'color', new_color);
                                
                    % Set line width and marker size
                    set(h_indiv(k), 'linewidth', 4, 'markersize', 15);
                    
                    % If condition name hasn't been added to legend yet
    %                 if ~any(strcmp(legend_names, data(i).condition{2}))
    %                     % Save the legend name based on the condition
    %                     legend_names{index} = data(i).condition{2};
    %                     h3(index) = h(i);
    %                     index = index + 1;
    %                 end

                    % Set properties of plot
                    set_plot_properties(h_indiv, ax_indiv, [], [], x_label,...% avg(i).condition{2}
                        y_label{j}, varargin{:});

                    if ~any(isnan(avg(i).blue))
                        % If plotting individual lines
                        if ~isempty(varargin{6}) && any(strcmp(varargin{6}, 'indiv'))
                            if any(strcmp(varargin{6}, 'indiv_separate'))
                                hold(ax_indiv, 'on');
                            end
                            close_to_line = false;
                            if close_to_line == true
                                temp_space = [10, 10, -10];
                                temp_y_blue = min(y(1:110),[],'all') - temp_space(k);
                            else
                                temp_y_blue = y_blue(i,j);
                            end

                            plot_blue_indiv = true;

                            if plot_blue_indiv == true
                                h_blue_indiv = plot(ax_indiv, [avg(i).blue(1,k),...
                                           avg(i).blue(2,k)],[temp_y_blue, temp_y_blue],...
                                           '.-', 'Color', [0 0.53 0.93, 0.5],...
                                           'MarkerEdgeColor', colors(c,:),...
                                           'Linewidth', 12,...
                                           'MarkerSize', 40);
                                uistack(h_blue_indiv,'bottom');
                            else
                                h_blue_avg = plot(ax_avg, [avg(i).avg_blue(1,b),...
                                   avg(i).avg_blue(2,b)],[y_blue(i,j), y_blue(i,j)],...
                                   '.-', 'Color', [0 0.53 0.93, 0.5],...
                                   'MarkerEdgeColor', colors(c,:),...
                                   'Linewidth', 12,...
                                   'MarkerSize', 40);
                                uistack(h_blue_avg,'bottom');
                            end
    
                            if any(strcmp(varargin{6}, 'indiv_separate'))
                                hold(ax_indiv, 'off');
                            end
                        end
                    end

                    h_indiv_legend(i) = h_indiv(k);
                end
                
                if ~isempty(varargin{6}) && any(strcmp(varargin{6}, 'indiv_separate'))
                    hold(ax_indiv, 'off');
                end
            end
            
            % If blue light was used
            if ~any(isnan(avg(i).avg_blue))
                % Plot blue bars
                for b = 1:size(avg(i).avg_blue, 2)
                    % If plotting average
                    if isempty(varargin{6}) || any(strcmp(varargin{6}, 'avg'))
                        h_blue_avg = plot(ax_avg, [avg(i).avg_blue(1,b),...
                                   avg(i).avg_blue(2,b)],[y_blue(i,j), y_blue(i,j)],...
                                   '.-', 'Color', [0 0.53 0.93, 0.5],...
                                   'MarkerEdgeColor', colors(c,:),...
                                   'Linewidth', 12,...
                                   'MarkerSize', 40);
                        uistack(h_blue_avg,'bottom');
                    end
                    
%                     % If plotting individual lines
%                     if ~isempty(varargin{6}) && any(strcmp(varargin{6}, 'indiv'))
%                         if any(strcmp(varargin{6}, 'indiv_separate'))
%                             hold(ax_indiv, 'on');
%                         end
% 
%                         h_blue_indiv = plot(ax_indiv, [avg(i).avg_blue(1,b),...
%                                    avg(i).avg_blue(2,b)],[y_blue(2,j), y_blue(2,j)],...
%                                    '.-', 'Color', [0 0.53 0.93, 0.5],...
%                                    'MarkerEdgeColor', colors(c,:),...
%                                    'Linewidth', 12,...
%                                    'MarkerSize', 40);
%                         uistack(h_blue_indiv,'bottom');
% 
%                         if any(strcmp(varargin{6}, 'indiv_separate'))
%                             hold(ax_indiv, 'off');
%                         end
%                     end
                end
            end
        end
        
        % If plotting average
        if isempty(varargin{6}) || (any(strcmp(varargin{6}, 'avg')) || any(strcmp(varargin{6}, 'indiv')))
            % Set properties of plot
            if isempty(varargin{6}) || any(strcmp(varargin{6}, 'avg'))
                set_plot_properties(h_avg, ax_avg, legend_names, [], x_label,...
                    y_label{j}, varargin{:});
            else
                set_plot_properties(h_indiv_legend, ax_indiv, legend_names, [], x_label,...
                        y_label{j}, varargin{:});
            end

            hold(ax_avg, 'off');
        else
        end
    end
end

function [colors, x_label, y_label, y_blue, field] = initiate_plot(n,...
                                                                  varargin)
%INITIATE_PLOT Initiates the plot
% 
%   Input
%       n: number of conditions
%       varargin: same varargin as entered into main function
%           5): enter 'norm' to normalize data by the first time point, or
%               'unnorm' to plot unnormalized data. The default is to plot
%               the normalized data.
%           6): enter 'indiv' to specify plotting individual lines. The
%               default is to plot the average
% 
%   Output
%       colors: array of colors to use
%       title_names: names for plots
%       field: field names for easily plotting average or individual data
%       y: y axis position of blue bars
% 
%   Overview
%       Defines the colors and picks the location of the blue light bars.
%       Change the value in y to manually change position of blue bars.
%       Also determines how plots should be made for average and individual
%       data.
    
    % Array of colors
    colors = [     0, 0.4470, 0.7410;   %blue
              0.8500, 0.3250, 0.0980;   %red
              0.9290, 0.6940, 0.1250;   %yellow/gold
              0.4940, 0.1840, 0.5560;   %purple
              0.4660, 0.6740, 0.1880;   %green
              0.3010, 0.7450, 0.9330;   %cyan
              0.6350, 0.0780, 0.1840;   %dark red
                   0,      0,      0;   %black
              50/256,71/256,200/256];

    
    
%     title_names = {'Mean # of Spots'};
%     % Titles for plots
%     title_names = {'Number of Spots',...
%                    'Average Area',...
%                    'Average Mean Intensity',...
%                    'Average Max Intensity',...
%                    'Average Sum of Intensities'};

    x_label = 'Time (min)';
    
    % If data isn't normalized
    if nargin >= 5 && strcmp(varargin{5}, 'unnorm')
        % If entered, set y-axis limits
        if nargin >= 8 && ~isempty(varargin{3})
            upper_y = varargin{3}(2);
        else
            upper_y = 200;
        end

        % y values close to x-axis at 0
        y_blue = [(upper_y / -25) .* (((upper_y / 400) + (upper_y / 360) * (n - 1)):-(upper_y / 360):(upper_y / 400)); % change 1 to 20
             (0.5 + 0.75 * (n - 1)):-0.75:0.5;
             325 .* ((1 + 0.6 * (n - 1)):-0.6:1);%500 .* ((0.5 + 0.75 * (n - 1)):-0.75:0.5);    2500 .* ((1.05 + 0.07 * (n - 1)):-0.07:1.05);
             500 .* ((0.5 + 0.75 * (n - 1)):-0.75:0.5);
             500 .* ((0.5 + 0.75 * (n - 1)):-0.75:0.5)]';

        y_label = {'# of Spots',...
               'Area (\mum^2)',...
               'Intensity (arb. unit)',...
               'Intensity (arb. unit)',...
               'Intensity (arb. unit)'};
    else
        % y values close to x-axis at 0
        y_blue = [-2.5 .* ((0.04 + 0.04 * (n - 1)):-0.04:0.04);
             0.6 + ((0.02 + 0.04 * (n - 1)):-0.04:0.02);
             0.7 + ((0.02 + 0.04 * (n - 1)):-0.04:0.02); %0.6 + ((0.02 + 0.04 * (n - 1)):-0.04:0.02);
             0.6 + ((0.02 + 0.04 * (n - 1)):-0.04:0.02);
             0.6 + ((0.02 + 0.04 * (n - 1)):-0.04:0.02)]';

        y_label = {'Number of Spots',...
               'Normalized Area (arb. unit)',...
               'Normalized Intensity (arb. unit)',...
               'Normalized Intensity (arb. unit)',...
               'Normalized Intensity (arb. unit)'};
    end

    field = {'n', 'A', 'avg_I', 'sum_I', 'max_I'};
%     field = {'avg_I'};
%     field = {'sum_I'};
%     field = {'A'};
%     field = {'max_I'};
%              'A';
%              'avg_I';
%              'max_I';
%              'sum_I'};
end

function set_plot_properties(h, ax, legend_names, title_names, x_label,...
    y_label, varargin)
%SET_PLOT_PROPERTIES Sets properties of the plots
%   
%   Inputs
%       h: handle to line objects for each plotted line object
%       legend_names: names of conditions
%       title_names: names of plots
%       varargin: same varargin as entered into main function
%           2): two entry array corresponding to the x-axis limits
%           3): two entry array corresponding to the y-axis limits
%       
%   Output
%       None
% 
%   Overview
%       Sets properties of the plots, such as legend, font size, and x and
%       y limits. x and y limits are determined from the inputs to the
%       function. Legend names are determined by the condition. The
%       location of the legend can be manually changed here.
    
    % If legend names are given
    if ~isempty(legend_names)
        % Make legend and set font size
        legend(h, legend_names, 'Location', 'northwest', 'Fontsize', 24)
        legend(ax, 'boxoff');
    end

    set(ax, 'Fontsize', 24);
    set(ax,'fontname','arial')

    if ~isempty(title_names)
        title(ax, title_names, 'FontSize', 24);
    end
    xlabel(ax, x_label, 'FontSize', 24);
    ylabel(ax, y_label, 'FontSize', 24);
    
%     % If normalized
%     if nargin >= 8 && ~strcmp(varargin{5}, 'unnorm')
%         set(gca, 'YScale', 'log');
%     end
    
    % If entered, set x-axis limits
    if nargin >= 7 && ~isempty(varargin{2})
        set(ax, 'xlim', varargin{2});
        set(ax,'XMinorTick','on')
    end
    
    % If entered, set y-axis limits
    if nargin >= 8 && ~isempty(varargin{3})
        set(ax, 'ylim', varargin{3});
    end
    
    ax.TickLength = [0.02, 0.02];
    ax.LineWidth = 2;
end

% function plot_boxplots(avg, varargin)
% %PLOT_BOX_PLOTS Plots boxplots
% %   
% %   Inputs
% %       h: handle to line objects for each plotted line object
% %       legend_names: names of conditions
% %       title_names: names of plots
% %       varargin: same varargin as entered into main function
% %           2): two entry array corresponding to the x-axis limits
% %           3): two entry array corresponding to the y-axis limits
% %       
% %   Output
% %       None
% % 
% %   Overview
% %       Sets properties of the plots, such as legend, font size, and x and
% %       y limits. x and y limits are determined from the inputs to the
% %       function. Legend names are determined by the condition. The
% %       location of the legend can be manually changed here.
% 
%     % Array of colors
%     colors = [     0, 0.4470, 0.7410;   %blue
%               0.8500, 0.3250, 0.0980;   %red
%               0.9290, 0.6940, 0.1250;   %yellow/gold
%               0.4940, 0.1840, 0.5560;   %purple
%               0.4660, 0.6740, 0.1880;   %green
%               0.3010, 0.7450, 0.9330;   %cyan
%               0.6350, 0.0780, 0.1840;   %dark red
%                    0,      0,      0];  %black
%     % Field names
%     field = {'max_n'; 'avg_n'; 'min_n'; 'sum_n'};
% 
%     % y labels
%     y_labels = {'Max # of Spots'; 'Mean # of Spots'; 'Min # of spots'; 'Area under the curve for # of Spots'};
% 
%     % x tick labels
%     x_tick_labels = {'nc12'; 'nc13'; 'nc14'};
% 
%     % For each field
%     for j = 1:size(field, 1)
%         % Initialize variables
%         size_arr = zeros(size(avg,2),1);
% 
%         % For each condition
%         for i = 1:size(avg, 2)
%             % Save the size of each replicate
%             size_arr(i) = size(avg(i).indiv,1);
%         end
%         
%         % Initialize variables for data, group names, and colors
%         y = NaN(max(size_arr), 3.* size(avg, 2));
%         g = string(3.* size(avg, 2));
%         c = zeros(3.* size(avg, 2), 3);
%         legend_names = cell(1, size(avg, 2));
% 
%         % For each condition
%         for i = 1:size(avg, 2)
%             condition_range = i:size(avg, 2):(i+2*size(avg, 2));
%             % Save 
%             y(1:size(avg(i).indiv,1), condition_range) = avg(i).(field{j}){:,:}; %(3*(i-1)+1):3*i
%             g(condition_range) = string(avg(i).(field{j}).Properties.VariableNames) + " " + avg(i).condition(2);
% 
%             % If color scheme entered is 2
%             if (nargin >= 5) && ~isempty(varargin{4}) && (varargin{4} >= 2)
%                 % Set color to color scheme to inputted number
%                 temp_c = avg(i).condition{2+varargin{4}};
%             else
%                 % Otherwise set color to color scheme 1
%                 temp_c = avg(i).condition{3};
%             end
%             
%             % Colors are replicated
%             [c(condition_range,:)] = repmat(colors(temp_c,:), [3,1]);
% 
%             % Save condition name for legend
%             legend_names{i} = avg(i).condition{2};
%         end
%         
%         % Initialize variables
%         x = zeros(size(y));
%         
%         % For each 
%         for k = 1:size(y, 2)
%             x(:,k) = k;
%         end
% 
%         % Make figure
%         figure;
% 
%         % Make a box plot of data y for groups g
%         boxplot(y, g, 'width', .75, 'colors', [0 0 0],...
%             'symbol', '');
% 
%         % Set formating for labels
%         set(gca, 'TickLabelInterpreter', 'tex');
% 
%         % Allow plotting without overwriting
%         hold on;
% 
%         % Plot scatter plot of real data
%         scatter(x, y, 100, c, 'filled',...
%                 'jitter', 'on', 'jitteramount', 0.2,...
%                 'MarkerFaceAlpha', .75,'MarkerEdgeAlpha', .75);
% 
%         % Turn off plotting without overwriting
%         hold off;
%         
%         % Set x tick locations
%         set(gca,'XTick', (size(avg, 2)-1):size(avg, 2):size(y, 2));
% 
%         % Label x ticks
%         xticklabels(x_tick_labels);
% 
%         % Label y
%         ylabel(y_labels{j}, 'FontSize', 18);
% 
%         % Make a legend
%         legend(legend_names, 'Location', 'northwest', 'Fontsize', 12);
%         legend('boxoff');
%     end
% end

function [stat, p] = plot_domain(avg, varargin)
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

    % Field names
%     field = {'sum_n'};
%     field = {'sum_n'; 'hull_area'};
    field = {'hull_area'};

    % y labels
%     y_labels = {'Area under the curve for # of Spots'; 'Area \mu^2'};

%     % x tick labels
%     x_tick_labels = {'nc12'; 'nc13'; 'nc14'};
    x_tick_labels = {'{\itdl-LEXY}'; '{\itdl-BLID}'};

    % Number of nuclear cycles to compare, which nc
    n = 1;
    nc = 3;

    stat = cell(size(field,2));
    p = cell(size(field,2));

    % For each field
    for j = 1:size(field, 1)

        % Initialize variables
        size_arr = zeros(size(avg,2),1);

        % For each condition
        for i = 1:size(avg, 2)
            % Save the size of each replicate
            size_arr(i) = size(avg(i).indiv,1);
        end
        
        % Initialize variables for data, group names, and colors
        y = NaN(max(size_arr), n.* size(avg, 2));
        g = strings(n.* size(avg, 2), 1);
        c = zeros(n.* size(avg, 2), 3);
        legend_names = cell(1, size(avg, 2));
        mean_y = zeros(1.* size(avg, 2), 1);
        err_y = zeros(1.* size(avg, 2), 1);

        % For each condition
        for i = 1:size(avg, 2)
            % Create a range for reordering data, based on the number of
            % nuclear cycles to plot
            condition_range = i:size(avg, 2):(i+(n-1)*size(avg, 2));
            
            % If data is to be normalized by subtraction
            if nargin >= 6 && strcmp(varargin{5}, 'norm')
                % Set the value to normalize by
                norm_val = avg(i).(field{j}){:,2};
            else
                norm_val = 0;
            end

            % Save the data as one array for individual points, mean and
            % error. Data is grouped by time, thus it is reordered to be
            % grouped by time first and then condition. The normalization
            % value is subtracted
            y(1:size(avg(i).indiv,1), condition_range) = avg(i).(field{j}){:,nc}-norm_val; %(3*(i-1)+1):3*i
            mean_y(condition_range) = mean(avg(i).(field{j}){:,nc}-norm_val, 1, 'omitnan');
            err_y(condition_range) = calc_error(avg(i).(field{j}){:,nc}-norm_val, 'SD', 1);
            
%             g(condition_range) = string(avg(i).(field{j}).Properties.VariableNames(nc)) + " " + avg(i).condition(2);
            g(condition_range) = avg(i).condition{1};

            % If color scheme entered is 2
            if (nargin >= 5) && ~isempty(varargin{4}) && (varargin{4} >= 2)
                % Set color of color scheme to inputted number
                temp_c = avg(i).condition{2+varargin{4}};
            else
                % Otherwise set color to color scheme 1
                temp_c = avg(i).condition{3};
            end
            
            % Colors are replicated
            [c(condition_range,:)] = repmat(colors(temp_c,:), [n,1]);

            % Save condition name for legend
            legend_names{i} = avg(i).condition{2};
        end
        
        [stat{j}, p{j}] = statistical_analysis(y, g);
        
        % Initialize variables
        x = zeros(size(y));
        
        % For each time point for each condition
        for k = 1:size(y, 2)
            x(:,k) = k;
        end

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
        
        % Plot the mean with the error bars and set properties
        h = errorbar(1:size(mean_y,1), mean_y, err_y, '.', 'Color', [0.63, 0.63, 0.63]);
        set(h, 'linewidth', 4, 'markersize', 50);

        % Turn off plotting without overwriting
        hold off;
        
        % Set properties of axis
        set(gca, 'ylim', [-7000, 8000],...
                 'ytick', -10000:2000:10000,...
                 'xlim', [0, size(mean_y,1)+1],...
                 'xtick', 1.5:2:size(mean_y,1),...%((size(avg, 2)+1)/2):size(avg, 2):size(y, 2),...
                 'fontsize', 30,...
                 'fontname', 'arial');

        % Label x ticks
%         xticklabels(x_tick_labels(nc));
        xticklabels(x_tick_labels);

        % Label y
%         ylabel(y_labels{j}, 'FontSize', 18);

        % Make a legend
        [~, hobj, ~, ~] = legend(unique(legend_names), 'Location', 'southwest', 'Fontsize', 30, 'fontname', 'arial');
        legend('boxoff');
        M = findobj(hobj,'type','patch');
        set(M,'MarkerSize',sqrt(500));
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