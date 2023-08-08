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

%     % Make a structure for storing data
%     data = struct('raw_image', cell(1,size(img_in,2)), 'name', [],...
%                   'nuc_mask', [], 'nuc_center', [],...
%                   'nuc_indices', [], 'embryo_center', [], 'tracks', [],...
%                   'areas', [], 'intensities', [], 'positions', [],...
%                   'spot_indices', []);

    % For each data set
    for i = 1:n
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
            data(i).spot_indices] = segment_nuclei(im(i).raw_image, im(i).threshold, data(i).t_align);
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

    reader = bfGetReader(path);
    omeMeta = reader.getMetadataStore();

    X = omeMeta.getPixelsSizeX(0).getValue();
    Y = omeMeta.getPixelsSizeY(0).getValue();
    Z = omeMeta.getPixelsSizeZ(0).getValue();
    T = omeMeta.getPixelsSizeT(0).getValue();
    C = omeMeta.getPixelsSizeC(0).getValue();
    
    % Allocate looped variable
    I = uint16(zeros(X,Y,C,Z));
    I2 = uint16(zeros(X,Y,C,T));
    raw_t = zeros(reader.getImageCount(),1);

    for t = 1:T
        for z = 1:Z
            for c = 1:C
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
        I2(:,:,:,t) = max(I, [], 4);
    end

    reader.close()
    
%     % Opens images using Bio-Formats for MATLAB
%     % https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html
%     im = bfopen(path);
    
%     % Save sizes of images in all dimensions, including time and color
%     % channels
%     X = im{1,4}.getPixelsSizeX(0).getValue();
%     Y = im{1,4}.getPixelsSizeY(0).getValue();
%     Z = im{1,4}.getPixelsSizeZ(0).getValue();
%     T = im{1,4}.getPixelsSizeT(0).getValue();
%     C = im{1,4}.getPixelsSizeC(0).getValue();
    
    % The physical length of a pixel
    xy_len = omeMeta.getPixelsPhysicalSizeX(0).value(...
                            ome.units.UNITS.MICROMETER); % in µm
    % Convert pixel length to a double
    xy_len = xy_len.doubleValue();

    % The physical distance between z slices
    z_len = omeMeta.getPixelsPhysicalSizeZ(0).value(...
                                ome.units.UNITS.MICROMETER); % in µm
    % Convert z length to a double
    z_len = z_len.doubleValue();
    
    % If image isn't projected, concatenate pixel and z length
    if strcmp(dims, '3D')
        vox_len = cat(2, xy_len, xy_len, z_len);
    % Otherwise, if it is projected, only concatenate pixel length for xy
    elseif strcmp(dims, '2D')
        vox_len = cat(2, xy_len, xy_len);
    end
    
%     % Allocate looped variable
%     raw_t = zeros(size(im{1,1},1),1);
    
%     % For each image plane
%     for i = 1:size(im{1,1},1)
%         % Try to get the time that elapsed during image aquisition
%         % If unable to, then the last z stack is incomplete
%         try
%             % Get time bewteen each z slice
%             raw_t(i,1) = im{1,4}.getPlaneDeltaT(0,...
%                                     i-1).value.doubleValue./60;
%             % Set to false since getting the time was a success
%             delete_last_t = false;
%         catch
%             % Set to true since getting the time was a failure
%             delete_last_t = true;
%         end
%     end
    
    % Reshape time to match the dimensions of channel, z, and time
    raw_t = reshape(raw_t,C,Z,T);
    
    % Times indicate when image aquisition finished, add preceding 0 to get
    % start of each z-stack and thus each time point
    t = [0;squeeze(raw_t(end,end,1:end-1))];

    % Reshape image data to match dimensions, X, Y, z, time, channels
    im = permute(I2, [1,2,5,4,3]);
    
%     % Reshape image data to match dimensions, X, Y, z, time, channels
%     im = permute(reshape(cat(3, im{1,1}{:,1}), Y, X, C, Z, T),...
%                   [1,2,4,5,3]);
    % Code to open large images in smaller parts
%     img1 = permute(reshape(cat(3, im{1,1}{1:(250.*C*Z),1}), Y, X, C, Z, 250),...
%                   [1,2,4,5,3]);
%     img1 = max(img1, [], 3);
%     img2 = permute(reshape(cat(3, im{1,1}{((250.*C*Z)+1):end,1}), Y, X, C, Z, T-250),...
%                   [1,2,4,5,3]);
%     img2 = max(img2, [], 3);
% 
%     im = cat(4, img1,img2);

%     % If 2D is desired
%     if strcmp(dims, '2D')
%         % Make a maximum z-projection
%         im = max(im, [], 3);
%     end
    
    % Delete last timepoint if z-stack is incomplete
    if delete_last_t
        im = im(:,:,:,1:(end-1),:);
        t = t(1:(end-1),1);
    end
end

function [mask, avg_I, mean_avg_I, centers, ind] = segment_nuclei(im, T, t_align, varargin)
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
%     ac = zeros(size(im, 4), 2);
%     embryo_center = zeros(size(im, 4), 2);
    
    % For each time point
    for i = 1:size(im,4)
        
%         % Gaussian blur
%         B = imgaussfilt(im(:,:,1,i,1), 20);
% 
%         % Threshold for entire embryo
%         bw = imbinarize(B, 0.05); %0.05
% 
%         % Morophologically close image to make a single object
%         se = strel('disk', 10);
%         J = imclose(bw, se);
% 
%         J = bwareaopen(J, 50000);
%         
%         % Get center of embryo
%         props = regionprops(J, 'Centroid');
% 
%         % To approximate mid line
%         % Take the difference of the mask, this will find edges
% %         dJ = diff(J, 1, 1);
%         dJ = diff(J, 1, 2);
%         q = sum(dJ, 2);
%         temp_i = false(512,511);
%         temp_i(:,end) = (q == 1);
%         dJ(temp_i) = -1;
%         temp_j = false(512,511);
%         temp_j(:,1) = (q == -1);
%         dJ(temp_j) = 1;
%         
%         if i == 177
%             dJ(5:6,508:511) = 0;
%         end
% 
%         % Get the row indices for the edges
% %         [r,~] = find(dJ);
%         [~,c] = find(dJ);
% 
%         % Reshape the row indices so each the first edge is one row and the
%         % second edge is the second row, for each pixel in the x dimension
% %         rr = reshape(r,2,size(im,2));
% 
%         cc = reshape(c,size(im,2),2);
% 
%         % Find the mid point of each column using average 
% %         mid = mean(rr,1);
%         mid = mean(cc,2);
% 
%         % Find the slopes using each point. Note the denominator is always
%         % one so it is left off (this is true because each column is an
%         % inreger and the difference is 1 between nearest columns)
%         m = diff(mid);
% 
%         % Find the average slope (to account for irregularites in the
%         % edges)
% %         avg_m = mean(m);
%         avg_m = 1./mean(m);
% 
%         % Save the slope for calculating the distance to the mid line
%         ac(i,1) = -avg_m;
% 
%         % Save the y-intercept for calculating the distance to the mid line
%         ac(i,2) = -(props.Centroid(2) - avg_m .* props.Centroid(1));
% 
%         % Save the center of the embryo
%         embryo_center(i,:) = props.Centroid;
        
        % Background subtract to remove any uneven illumination
%         im_bg_subtract = im(:,:,1,i,1) - imgaussfilt(im(:,:,1,i,1), 50);
%         im_bg_subtract = im(:,:,1,i,1) - imgaussfilt(im(:,:,1,i,1), 50);

        % Blur for segmenting nuclei
        B = imgaussfilt(im(:,:,1,i,1), 4);

        % Threshold image for nuclei, note that threshold is different
        % because of background subtraction
%         bw = imbinarize(B, .1);

        if isempty(T)
            T = 0.1;
        end

        if isempty(t_align)
            t_align = 100;
        end

        bw = imbinarize(B, T+((i-1) * (.00125/t_align)));

        % Morophologically open image to disconnect nuclei
        se = strel('disk', 5);
        J = imopen(bw, se);

        se = strel('disk', 3);
        J = imclose(J, se);
        J = imfill(J,'holes');

        % Remove objects/nuclei touching the edge of the image
        J = imclearborder(J);

        D = bwdist(~J);
        J = imhmin(-D,1);
        L = watershed(J);
        L(~bw) = 0;
        L1 = logical(L);
        mask(:,:,1,i) = bwareafilt(L1,[100,700]);

        % Remove small and large objects
%         mask(:,:,1,i) = bwareafilt(J, [300,1200]);
%         mask(:,:,1,i) = bwareafilt(J, [100,300]);

        % Watershed
%         B = imgaussfilt(im(:,:,1,i), 5);
%         bw = imbinarize(B, .045);
%         D = bwdist(~bw);
%         J = imhmin(-D,1);
%         L = watershed(J);
%         L(~bw) = 0;
%         L1 = logical(L);
%         mask(:,:,1,i) = bwareafilt(bw2,[300,1200]);

        % Get properties for nuclei including the center, and a list of
        % pixels in each object
        props = regionprops(mask(:,:,1,i), im(:,:,1,i), 'Centroid',...
            'MeanIntensity', 'PixelIdxList');
        
        avg_I{i} = cat(1, props.MeanIntensity);
        
        if ~isempty(avg_I{i})
            mean_avg_I(i,1) = mean(avg_I{i}, 1);
            mean_avg_I(i,2) = calc_error(avg_I{i}, err_method, 1);
        end

        % Save centers of nuclei
        centers{i} = cat(1, props.Centroid);

        % Save list of pixels
        ind{i} = {props.PixelIdxList}';

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