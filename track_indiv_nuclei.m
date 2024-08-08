function [track_data] = track_indiv_nuclei(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    % Make a structure for storing data
    track_data = struct('name', cell(1,size(data,2)), 'time', [],...
                  'centers', [], 'I', [], 't_blue_on', [],...
                  't_blue_off', [], 'params_on', [], 'params_off', [],...
                  'mean_params_on', [], 'mean_params_off', []);

    % For each data set
    for i = 1:size(data,2)
        % Save name
        track_data(i).name = data(i).name;
        track_data(i).time = data(i).time;

        % Find and save index of closest cells over time a certain
        % distance from the mid line, defined by constants in ac
        [track_data(i).centers, track_data(i).I] = index_of_closest_cell(data(i).centers, data(i).avg_I);
    end
end

function [tracks_centers, tracks_I] = index_of_closest_cell(centers, I)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     % Initialize variables
%     sizes = zeros(size(nuc_center,1),1);
%     
%     % For each time point, find the number of spots
%     for n = 1:size(nuc_center,1)
%         sizes(n) = size(nuc_center{n,1},1);
%     end
    
    % Intialize tracks. 100 is an overestimate of number of nuclei. Other
    % dimension is the number of time points. The third dimension is for
    % the centers and the intensity.
    tracks_centers = nan(100, 2, size(centers,1));
    tracks_I = nan(100, size(centers,1));
    n_tracks = size(centers{1,1},1);
    tracks_centers(1:n_tracks,:,1) = centers{1,1};
    tracks_I(1:n_tracks,1) = I{1,1};
    
    % For each time point
    for t = 2:size(centers,1)
        % For each nuclei
        for n = 1:size(centers{t,1},1)
            for i = t-1:-1:1
                if i < t-10
                    n_tracks = n_tracks + 1;
                    tracks_centers(n_tracks,:,t) = centers{t,1}(n,:);
                    tracks_I(n_tracks,t) = I{t,1}(n,1);
                    break
                end
                d = calc_dist(tracks_centers(:,:,i), centers{t,1}(n,:));
    
                % Save the index of the shortest distance
                [min_d, ind] = min(d);
    
                if min_d < 25
                    tracks_centers(ind,:,t) = centers{t,1}(n,:);
                    tracks_I(ind,t) = I{t,1}(n,1);
                    break
                elseif i == 1
                    n_tracks = n_tracks + 1;
                    tracks_centers(n_tracks,:,t) = centers{t,1}(n,:);
                    tracks_I(n_tracks,t) = I{t,1}(n,1);
                end
            end
        end
    end

    for m = size(tracks_I,1):-1:1
        if (sum(isnan(tracks_I(m,:))) > size(centers,1) / 2) || (all(isnan(tracks_I(m,end-10:end)))) || (all(isnan(tracks_I(m,1:10))))
            tracks_centers(m,:,:) = [];
            tracks_I(m,:) = [];
        end
    end
end

function d = calc_dist(x, y)
%CALC_DIST Calculates the distance between two points in n-dimensions
%   
%   Input
%       x: corrdinates for an array of points (can be mxn in size)
%       y: corrdinates for a point (should be 1xn in size)
%   
%   Output
%       d: array of distances bewteen the points in x and the point in y
%
%   Overview
%       Calculates the distance in n-dimensial space by subtracting x and y
%       by applying element-wise operation to the two arrays with implicit
%       expansion enabled. This will subtract y from each row of x if x is
%       mxn and y is 1xn, where n is the number of dimensions. These values
%       are then squared and summed upon the second dimension of the array.
%       Finally the square root is taken. This gives the distance formula,
%       d = sqrt((x1-x2)^2+(y1-y2)^2) but for n-dimensions and for
%       mutiple points in x from one point y. Note that if x and y are both
%       mxn, then each row in y will be subtracted from the corresponding
%       row in x.

    d = sqrt(sum((x-y).^2, 2));
end