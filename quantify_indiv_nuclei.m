function [im, data] = quantify_indiv_nuclei(varargin)
%UNTITLED7 Summary of this function goes here
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

        if size(varargin,1) >=2 && ~isempty(varargin{2})
            specify_z = varargin{2};
        else
            specify_z = [];
        end

        im = struct('folder', cell(1,n),...
                    'name', [],...
                    'condition', [],...
                    'signal_channels', 1,...
                    'background_channel', [],...
                    'raw_image', [],...
                    'threshold', [],...
                    'mask', [],...
                    'image_dims', []);

        data = struct('folder', cell(1,n),...
                    'name', [],...
                    'pixel_length', [],...
                    'time', [],...
                    'raw_time', [],...
                    'avg_I', [],...
                    'centers', NaN,...
                    'spot_indices', NaN,...
                    'ind', NaN,...
                    'blue_light', NaN,...
                    't_align', 1,...
                    't_norm', NaN,...
                    'params', []);

        name = cell(n,1);
        folder = cell(n,1);
        for i = 1:n
            % Use menu to select file
            [name{i},folder{i}] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
                'Select the microscope images', 'MultiSelect', 'off');
        end
    end

    % For each data set
    for i = 1:n
        % If a structure was not the first input
        if ~isstruct(varargin{1})
            % Open image/movie
            [im(i).folder, im(i).name, im(i).raw_image,...
                data(i).pixel_length, data(i).time,...
                data(i).raw_time] = open_img('2D', name{i}, folder{i}, specify_z);
    
            % Copy folder and filename to data structure
            data(i).folder = im(i).folder;
            data(i).name = im(i).name;
        end

        data(i).avg_I = cell(size(im(i).raw_image,4),1);
        data(i).centers = cell(size(im(i).raw_image,4),1);

        im(i).mask = false(size(im(i).raw_image(:,:,:,:,1)));
    
        for t = 1:size(im(i).raw_image,4)
            im(i).mask(:,:,1,t) = segment_nuclei(im(i).raw_image(:,:,:,t,1));

            if t == 13 || ((t >= 40) && (t < 47))
                im(i).mask(:,:,1,t) = im(i).mask(:,:,1,t-1);
            end
            
            props = regionprops(im(i).mask(:,:,1,t), im(i).raw_image(:,:,:,t,1), 'Centroid', 'MeanIntensity','PixelIdxList');
            data(i).avg_I{t} = cat(1, props.MeanIntensity);
            data(i).centers{t} = cat(1, props.Centroid);
        end
    end
end

function [path, embryo_number, im, vox_len, t, raw_t] = open_img(dims, name, folder, specify_z)
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
    
%     % Use menu to select file
%     [name,folder] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
%         'Select the microscope images', 'MultiSelect', 'off');
    
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

        if isempty(specify_z) || strcmp(specify_z, 'proj')
            % Make a max intensity projection
%             I2(:,:,:,t) = max(I, [], 4);
            I2(:,:,:,t) = mean(I,4);
        elseif isfloat(specify_z) && ~isempty(specify_z)
            I2(:,:,:,t) = I(:,:,:,specify_z);
        end
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

function bw = segment_nuclei(im)
%AXIS_POINTS Determine points along major axis of an ellipse.
%   
%   Input
%       img: the raw image, a z-projection, returned from open_img
%   
%   Output
%       p: positions for centering ROIs along major axis of an ellipse,
%           determined from the shape of a drosophila embryo
%
%   Overview
%       Segments an image of a drosophila embryo, using Otsu's method for
%       thresholding, then gets the ellipse of the mask. The function uses
%       the properties of the ellipse to calculate different positions
%       along the AP (major axis) of the embryo. ROIs are made with these
%       points at the center. A line perpendicular to the axis is drawn
%       in case the ROI has to be moved due to rotation of the embryo.
%       Positions include: 10%, 20%, 25%, 30%, 40%, 50%, 60%, 70%, 75%,
%       80%, 90%. Size of the ROI may make these positions mean less if the
%       ROI is large, either due to inclusion of background pixels, or
%       highly overlapping domains.

    im_LoG = edge(im,'log',0,10);
    im_LoG = imfill(im_LoG, 'holes');
    se = strel('disk', 5);
    bw = imopen(im_LoG, se);

    bw = bwareafilt(bw, [1000, 10000]);

    do_watershedding = true;

    if do_watershedding
        % For watershedding, find the distances in the mask
        D = bwdist(~bw);

        % Only keep certain minimums
        J = imhmin(-D,1);

        % Perform watershed
        L = watershed(J);

        % Remove mask of wathershed that's outside of the original mask
        L(~bw) = 0;

        % Convert to logical and filter out small and large objects
        bw = logical(L);

        bw = bwareafilt(bw, [1000, 10000]);
    end
end