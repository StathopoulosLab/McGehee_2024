function [t_out, raw_t_out, name] = get_time_from_image(n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% new_t = t-t(92);
% D = duration(minutes(new_t),'format','hh:mm');
    
    
    t_out = cell(n,1);
    raw_t_out = cell(n,1);
    name = cell(n,1);
    folder = cell(n,1);

    for j = 1:n
        % Use menu to select file
        [name{j}, folder{j}] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
            'Select the microscope images', 'MultiSelect', 'off');
    end

    for k = 1:n
        % Construct full path
        path = fullfile(folder{k}, name{k});
        
        % Opens images using Bio-Formats for MATLAB
        % https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html
    %     img = bfopen(path);
    
        % Construct a Bio-Formats reader
        reader = bfGetReader(path);
        omeMeta = reader.getMetadataStore();
        
        % Save sizes of images in all dimensions, including time and color
        % channels
        Z = omeMeta.getPixelsSizeZ(0).getValue();
        T = omeMeta.getPixelsSizeT(0).getValue();
        C = omeMeta.getPixelsSizeC(0).getValue();
        
        % Allocate looped variable
        raw_t = zeros(reader.getImageCount(),1);
    
        for i = 1:reader.getImageCount()
            % Try to get the time that elapsed during image aquisition
            % If unable to, then the last z stack is incomplete
            try
                % Get time bewteen each z slice
                raw_t(i,1) = omeMeta.getPlaneDeltaT(0,...
                                        i-1).value.doubleValue./60;
                % Set to false since getting the time was a success
                delete_last_t = false;
            catch
                % Set to true since getting the time was a failure
                delete_last_t = true;
            end
        end
        
        % Close the reader
        reader.close()
    
        % Reshape time to match the dimensions of channel, z, and time
        raw_t = reshape(raw_t,C,Z,T);
        
        % Times indicate when image aquisition finished, add preceding 0 to get
        % start of each z-stack and thus each time point
        t = [0;squeeze(raw_t(end,end,1:end-1))];
    
        % Delete last timepoint if z-stack is incomplete
        if delete_last_t
            t = t(1:(end-1),1);
        end

        t_out{k} = t;
        raw_t_out{k} = raw_t;
    end
end