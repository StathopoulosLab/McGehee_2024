function [im, data] = quantify_nuclei(varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    % If first input is a structure
    if isstruct(varargin{1})
        % n is size of structure, first input is the image data, second is 
        % the data structure
        n = size(varargin{1}, 2);
        im = varargin{1};
        data = varargin{2};
    else
        % n is first input and make the data structures for output
        n = varargin{1};

        im = struct('folder', cell(1,n),...
                    'name', [],...
                    'condition', [],...
                    'signal_channels', [],...
                    'background_channel', [],...
                    'raw_image', [],...
                    'threshold', [],...
                    'mask', [],...
                    'image_dims', []);
    
        data = struct('folder', cell(1,n),...
                    'name', [],...
                    'condition', [],...
                    'pixel_length', [],...
                    'time', [],...
                    'raw_time', [],...
                    'avg_I', [],...
                    'mean_avg_I', [],...
                    'centers', NaN,...
                    'spot_indices', NaN,...
                    'ind', NaN,...
                    'blue_light', NaN,...
                    't_align', 1,...
                    't_norm', NaN,...
                    'nuc_cycle', [NaN, NaN; NaN, NaN; NaN, NaN], ...
                    'pts', [],...
                    'hull_pts', [],...
                    'rm_pts', []);
    end

    % For each data set
    for i = 1:n
        % If a structure was not the first input
        if ~isstruct(varargin{1})
            % Open image/movie
            [im(i).folder, im(i).name, im(i).raw_image,...
                data(i).pixel_length, data(i).time,...
                data(i).raw_time] = open_img('2D');
    
            % Copy folder and filename to data structure
            data(i).folder = im(i).folder;
            data(i).name = im(i).name;
        end

        % Segment nuclei using background fluorescence
        [im(i).mask, data(i).avg_I, data(i).mean_avg_I, data(i).centers,...
            data(i).spot_indices] = segment_nuclei(im(i).raw_image,...
            im(i).threshold, data(i).t_align);
    end

end

function [path, embryo_number, im, vox_len, t, raw_t] = open_img(dims)
%OPEN_IMG Open a czi with a z-stack, a time series, and channels
% 
%   Inputs
%       dims: '2D' or '3D' to determine if a z-projection is made or not
% 
%   Outputs
%       path: the folder path containing the opened file
%       embryo_number: part of the file name before the first space
%       im: raw images or maximum z-projection of images
%       t: the time for each z-projection, taken as the time of the last
%           z-plane in a z-stack
%       raw_t: the time to finish each z-plane
% 
%   Overview
%       Opens the selected image file. The last time point is deleted if
%       the final z-stack is not complete.
    
    % Use menu to select file
    [name,folder] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
        'Select the microscope images', 'MultiSelect', 'off');
    
    % Construct full path
    path = fullfile(folder,name);
    
    % Split and save part of file name before first space as unique
    % identifier
    file_name_parts = strsplit(name, ' ');
    embryo_number = file_name_parts{1};
    
    % Use bioformats to read in file
    reader = bfGetReader(path);
    omeMeta = reader.getMetadataStore();

    % Save the size of X, Y, Z, T, and C
    X = omeMeta.getPixelsSizeX(0).getValue();
    Y = omeMeta.getPixelsSizeY(0).getValue();
    Z = omeMeta.getPixelsSizeZ(0).getValue();
    T = omeMeta.getPixelsSizeT(0).getValue();
    C = omeMeta.getPixelsSizeC(0).getValue();
    
    % Allocate looped variable
    I = uint16(zeros(X,Y,C,Z));
    I2 = uint16(zeros(X,Y,C,T));
    raw_t = zeros(reader.getImageCount(),1);
    
    % For each time step
    for t = 1:T
        % For each z slice
        for z = 1:Z
            % For each channel
            for c = 1:C
                % Get the index and save the image
                i = reader.getIndex(z-1, c-1, t-1)+1;
                I(:,:,c,z) = bfGetPlane(reader, i);

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
        end
        % Make a mex intensity projection
        I2(:,:,:,t) = max(I, [], 4);
    end
    
    % Close the open file
    reader.close()
    
    % Save the physical length of a pixel in µm
    xy_len = omeMeta.getPixelsPhysicalSizeX(0).value(...
                            ome.units.UNITS.MICROMETER);
    
    % Convert pixel length to a double
    xy_len = xy_len.doubleValue();

    % Save the physical distance between z slices in µm
    z_len = omeMeta.getPixelsPhysicalSizeZ(0).value(...
                                ome.units.UNITS.MICROMETER);

    % Convert z length to a double
    z_len = z_len.doubleValue();
    
    % If image isn't projected, concatenate pixel and z length
    if strcmp(dims, '3D')
        vox_len = cat(2, xy_len, xy_len, z_len);
    % Otherwise, if it is projected, only concatenate pixel length for xy
    elseif strcmp(dims, '2D')
        vox_len = cat(2, xy_len, xy_len);
    end
    
    % Reshape time to match the dimensions of channel, z, and time
    raw_t = reshape(raw_t,C,Z,T);
    
    % Times indicate when image aquisition finished, add preceding 0 to get
    % start of each z-stack and thus each time point
    t = [0;squeeze(raw_t(end,end,1:end-1))];

    % Reshape image data to match dimensions, X, Y, z, time, channels
    im = permute(I2, [1,2,5,4,3]);
    
    % Delete last timepoint if z-stack is incomplete
    if delete_last_t
        im = im(:,:,:,1:(end-1),:);
        t = t(1:(end-1),1);
    end
end

function [mask, avg_I, mean_avg_I, centers, ind] = segment_nuclei(im, T,...
    t_align, varargin)
%SEGMENT_NUCLEI Segments nuclei based on MCP-GFP/RFP background 
%   The function

    if isempty(varargin)
        err_method = 'SEM';
    else
        err_method = varargin{1};
    end
    
    % Initialize variables
    mask = false(size(im(:,:,:,:,1)));
    avg_I = cell(size(im, 4), 1);
    mean_avg_I = zeros(size(im, 4), 2);
    centers = cell(size(im, 4), 1);
    ind = cell(size(im, 4), 1);
    
    % For each time point
    for i = 1:size(im,4)
        % Blur for segmenting nuclei
        B = imgaussfilt(im(:,:,1,i,1), 4);
    
        % If no threshold is provided set it to 0.1
        if isempty(T)
            T = 0.1;
        end
        
        % If no t_align is provided set it to 100.
        if isempty(t_align)
            t_align = 100;
        end

        % Threshold image, increasing the threshold every image
        bw = imbinarize(B, T+((i-1) * (.00125/t_align)));

        % Morophologically open image to disconnect nuclei
        se = strel('disk', 5);
        J = imopen(bw, se);
        
        % Morphologically close image and fill to remove any holes
        se = strel('disk', 3);
        J = imclose(J, se);
        J = imfill(J,'holes');

        % Remove objects/nuclei touching the edge of the image
        J = imclearborder(J);
        
        % For wathershedding, find the distances in the mask
        D = bwdist(~J);

        % Only keep certain minimums
        J = imhmin(-D,1);

        % Perform watershed
        L = watershed(J);

        % Remove mask of wathershed that's outside of the original mask
        L(~bw) = 0;

        % Convert to logical and filter out small and large objects
        L1 = logical(L);
        mask(:,:,1,i) = bwareafilt(L1,[100,700]);

        % Get properties for nuclei including the center, a list of pixels
        % in each object, and mean intensity of object
        props = regionprops(mask(:,:,1,i), im(:,:,1,i), 'Centroid',...
            'MeanIntensity', 'PixelIdxList');
        
        % Save the mean intensities
        avg_I{i} = cat(1, props.MeanIntensity);
        
        % If there are mean intensities, calculate the mean of mean
        % intensities
        if ~isempty(avg_I{i})
            mean_avg_I(i,1) = mean(avg_I{i}, 1);
            mean_avg_I(i,2) = calc_error(avg_I{i}, err_method, 1);
        end

        % Save centers of nuclei
        centers{i} = cat(1, props.Centroid);

        % Save list of pixels
        ind{i} = {props.PixelIdxList}';
        
        % If less than 10 objects were detected it is unlikely those are
        % real nuclei, so remove all detected data
        if size(avg_I{i}) < 10
            avg_I{i} = [];
            centers{i} = [];
            mean_avg_I(i,:) = [0,0];
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