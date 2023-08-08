function [n_spots, xy, qn_spots, new_xy, path] = import_spot_analysis
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    % Open file
    [n_spots, xy, path] = open_excel;

    % Save only the number of spots
    qn_spots = n_spots{:,2};
    
    new_xy = cell(size(xy,1),1);

    for i = 1:size(xy,1)
        new_xy{i} = reshape(xy{i,2:end}, 2, [])';
    end
end

function [n_spots, xy, path] = open_excel
%OPEN_EXCEL Open a excel file wih Activations and Coords X-Y sheets
% 
%   Inputs
%       dims: '2D' or '3D' to determine if a z-projection is made or not
% 
%   Outputs
%       n_spots: number of spots at a time point
%       xy: xy coordinates
% 
%   Overview
%       Opens the selected excel file.
    
    % Use menu to select file
    [name,folder] = uigetfile({'*.xlsx', 'Excel files (*.xlsx)'},...
        'Select one file', 'MultiSelect', 'off');
    
    % Construct full path
    path = fullfile(folder,name);
    
    % Read in Activations sheet and Coords X-Y sheet
    n_spots = readtable(path, 'Sheet', 'Activations');
    xy = readtable(path, 'Sheet', 'Coords X-Y');
end