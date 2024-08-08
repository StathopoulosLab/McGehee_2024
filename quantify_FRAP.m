function [data] = quantify_FRAP(varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    % If first input is a structure
    if isstruct(varargin{1})
        % n is size of structure, first input is the image data, second is 
        % the data structure
        n = size(varargin{1}, 2);
        data = varargin{1};
    else
        % n is first input and make the data structures for output
        n = varargin{1};
    
        data = struct('folder', cell(1,n),...
                    'name', [],...
                    'raw_image', [],...
                    'roi_mask', [],...
                    'quant_mask', [],...
                    'mask', [],...
                    'single_mask', [],...
                    'pixel_length', [],...
                    'time', [],...
                    'raw_time', [],...
                    'avg_I', [],...
                    'roi_center', [],...
                    'roi_radii', zeros(1,2),...
                    'centers', NaN,...
                    'spot_indices', NaN,...
                    'ell_axes', NaN,...
                    'circularity', NaN,...
                    'ind', NaN,...
                    'blue_light', NaN,...
                    't_align', 1,...
                    't_norm', NaN,...
                    'nuc_cycle', [NaN, NaN; NaN, NaN; NaN, NaN], ...
                    'pts', [],...
                    'hull_pts', [],...
                    'rm_pts', [],...
                    't_bleach', [],...
                    'params', []);

        name = cell(n,1);
        folder = cell(n,1);
        for i = 1:n
            % Use menu to select file
            [name{i},folder{i}] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
                'Select the microscope images', 'MultiSelect', 'off');
        end
    end

    for i = 1:n
        if ~isstruct(varargin{1})
            % Open image/movie
            [data(i).folder, data(i).name, data(i).raw_image,...
                data(i).pixel_length, data(i).time, data(i).raw_time,...
                data(i).roi_center, data(i).roi_radii] = open_img(name{i}, folder{i});
        end
    
%         figure;
%         imshow(data(i).raw_image(:,:,1,1,1), [0 15000]);
    
%         h1 = drawellipse('Center',data(i).roi_center,'Semiaxes',data(i).roi_radii,'Color','g');
    
%         h2 = drawellipse('Color','y');
%         h2 = customWait(h2);


%         data(i).roi_mask = createMask(h1,data(i).raw_image(:,:,:,1,1));
%         data(i).quant_mask = createMask(h2,data(i).raw_image(:,:,:,1,1));
        
        data(i).avg_I = zeros(size(data(i).raw_image,4), 1);
        data(i).centers = zeros(size(data(i).raw_image,4), 2);
        
        data(i).mask = false(size(data(i).raw_image(:,:,:,:,1)));
        data(i).single_mask = false(size(data(i).raw_image(:,:,:,:,1)));
        centers = cell(size(data(i).raw_image,4),1);
    
        for t = 1:size(data(i).raw_image,4)
            data(i).mask(:,:,1,t) = segment_nuclei(data(i).raw_image(:,:,:,t,1));
            
%             if (t > 5) && (t < 15)
%                 data(i).mask(:,:,1,t) = data(i).mask(:,:,1,5);             
%             elseif (t > 20) && (t < 125)
%                 data(i).mask(:,:,1,t) = data(i).mask(:,:,1,20);
%             end

            props = regionprops(data(i).mask(:,:,1,t), data(i).raw_image(:,:,:,t,1), 'Centroid', 'MeanIntensity','PixelIdxList');
            centers{t} = cat(1, props.Centroid);
            
            if (t ~= 1) && (size(centers{t}, 1) < 30)
                data(i).mask(:,:,1,t) = data(i).mask(:,:,1,t-1);
                props = regionprops(data(i).mask(:,:,1,t), data(i).raw_image(:,:,:,t,1), 'Centroid', 'MeanIntensity','PixelIdxList');
                centers{t} = cat(1, props.Centroid);
            end
            
            if t == 1
                [min_d,min_ind] = min(calc_dist(centers{t}, data(i).roi_center));
            else
                [min_d,min_ind] = min(calc_dist(centers{t}, data(i).centers(t-1,:)));
            end

            if min_d < 25
                temp_mask = false(size(data(i).raw_image,1),size(data(i).raw_image,2));
                temp_mask(props(min_ind).PixelIdxList) = true;
                data(i).single_mask(:,:,1,t) = temp_mask;
                data(i).avg_I(t) = props(min_ind).MeanIntensity;
                data(i).centers(t,:) = props(min_ind).Centroid;
            elseif t ~= 1
                data(i).single_mask(:,:,1,t) = data(i).single_mask(:,:,1,t-1);
                prop_single = regionprops(data(i).single_mask(:,:,1,t), data(i).raw_image(:,:,:,t,1), 'Centroid', 'MeanIntensity');
                data(i).avg_I(t) = prop_single.MeanIntensity;
                data(i).centers(t,:) = prop_single.Centroid;
            end
        end

%         if i >=2
%             t{i} = t{i} + mean(diff(t{i}),1) + t{i-1}(end);
%         end
    end
end

function [path, file_name, im, pix_len, time, raw_t, center, r_xy] = open_img(name, folder)
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
    
    center = zeros(1,2);
    r_xy = zeros(1,2);
    

%     % Use menu to select file
%     [name,folder] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
%         'Select the microscope images', 'MultiSelect', 'off');
    
    % Construct full path
    path = fullfile(folder,name);

    % Split and save part of file name before first space as unique
    % identifier
    file_ext = strsplit(name, '.');

    file_name = file_ext{1};

    % Use bioformats to read in file
    reader = bfGetReader(path);
    omeMeta = reader.getMetadataStore();

    % Save the size of X, Y, Z, T, and C
    X = omeMeta.getPixelsSizeX(0).getValue();
    Y = omeMeta.getPixelsSizeY(0).getValue();
    Z = omeMeta.getPixelsSizeZ(0).getValue();
    T = omeMeta.getPixelsSizeT(0).getValue();
    C = omeMeta.getPixelsSizeC(0).getValue();
    
    r_xy(1,1) = double(omeMeta.getEllipseRadiusX(0,0));
    r_xy(1,2) = double(omeMeta.getEllipseRadiusY(0,0));
    center(1,1) = double(omeMeta.getEllipseX(0,0));
    center(1,2) = double(omeMeta.getEllipseY(0,0));
    
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
        % Make a max intensity projection
%         I2(:,:,:,t) = max(I, [], 4);
        I2(:,:,:,t) = mean(I,4);
    end
    
    % Close the open file
    reader.close()
    
    % Save the physical length of a pixel in µm
    xy_len = omeMeta.getPixelsPhysicalSizeX(0).value(...
                            ome.units.UNITS.MICROMETER);
    
    % Convert pixel length to a double
    xy_len = xy_len.doubleValue();

    pix_len = cat(2, xy_len, xy_len);
    
    % Reshape time to match the dimensions of channel, z, and time
    raw_t = reshape(raw_t,C,Z,T);
    
    % Times indicate when image aquisition finished, add preceding 0 to get
    % start of each z-stack and thus each time point
    time = [0;squeeze(raw_t(end,end,1:end-1))];

    % Reshape image data to match dimensions, X, Y, z, time, channels
    im = permute(I2, [1,2,5,4,3]);
    
    % Delete last timepoint if z-stack is incomplete
    if delete_last_t
        im = im(:,:,:,1:(end-1),:);
        time = time(1:(end-1),1);
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
%       mxn and y is 1xn, where n is the number of dimensions. These valued
%       are then squared and summed upon the second dimension of the array.
%       Finally the square root is taken. This gives the distance formula,
%       d = sqrt((x1-x2)^2+(y1-y2)^2) but for n-dimensions and for
%       mutiple points in x from one point y.

    d = sqrt(sum(bsxfun(@minus, x, y).^2,2));
end   