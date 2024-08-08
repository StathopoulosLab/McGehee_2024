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
                    'condition', [],...
                    'pixel_length', [],...
                    'time', [],...
                    'raw_time', [],...
                    'avg_I', [],...
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
                    'rm_pts', []);

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

        % Segment nuclei using background fluorescence
        [im(i).mask, data(i).avg_I, data(i).centers,...
            data(i).spot_indices, data(i).ell_axes, data(i).circularity] =...
            process_im(im(i).raw_image, im(i).signal_channels);
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

function [mask, avg_I, centers, ind, ellipse_axes, circularity] = ...
    process_im(im, ch)
%SEGMENT_NUCLEI Segments nuclei based on MCP-GFP/RFP background 
%   The function
    
    % Initialize variables
    mask = false(size(im(:,:,:,:,ch)));
    avg_I = cell(size(im, 4), 1);
    centers = cell(size(im, 4), 1);
    ellipse_axes = cell(size(im, 4), 1);
    ind = cell(size(im, 4), 1);
    circularity = cell(size(im, 4), 1);
%     avg_I_em = cell(size(im, 4), 1);
    
    % For each time point
    for t = 1:size(im,4)
        bw = segment_nuclei(im(:,:,1,t,ch));

        % Get properties for nuclei including the center, a list of pixels
        % in each object, and mean intensity of object
        props = regionprops(bw, im(:,:,1,t,ch), 'Centroid',...
            'PixelIdxList');
        
        % Save centers of nuclei
        centers{t} = cat(1, props.Centroid);

        % Save list of pixels
        ind{t} = {props.PixelIdxList}';

        decrease_thresh = true;
        T_ind = 1:(1/size(im,4)):3;

        if decrease_thresh
%             T = 0.02 / T_ind(t);
            T = 0.005;
        end

        [em_bw, boundary] = segment_embryo(im(:,:,1,t,ch), T);

        if ~isempty(boundary)
            size_im = size(im,1);
            af = fit_ellipse([size_im/2, size_im/2, 2*size_im, size_im, pi/4], boundary(:,2), boundary(:,1), size_im);
    
            p = zeros(2);
            % Calculate positions along right side of major axis
            p(1,:) = parameterized_ellipse(0, af(3), af(4), af(5), [af(1), af(2)]);
            % Calculate positions along left side of major axis
            p(2,:) = parameterized_ellipse(0, -af(3), af(4), af(5), [af(1), af(2)]);

            for i = 1:size(centers{t},1)
                if calc_dist_to_line(p, centers{t}(i,:)) > 100
                    bw(ind{t}{i}) = false;
                end
            end
        end

        

%         figure;
%         ax = axes;
%         hold on;
%         imshow(im(:,:,1,t,ch), [0, 10000], 'Parent', ax);
%         scatter(ax, boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
%         phi = 0:0.01:2*pi;
%         new_xy = parameterized_ellipse(phi, af(3), af(4), af(5), [af(1), af(2)]);
%         plot(ax, new_xy(:,1), new_xy(:,2), 'white', 'LineWidth', 3);
%         plot(ax, cat(1, p(1,1), af(1), p(2,1)), cat(1, p(1,2), af(2), p(2,2)), 'red');
%         scatter(ax, centers{t}(:,1), centers{t}(:,2), 'b');
%         hold off
%         pause;
        
        

        check_bound = true;

        if check_bound
            bw(~em_bw) = false;
        end

%         % Get properties for nuclei including the center, a list of pixels
%         % in each object, and mean intensity of object
%         props_em = regionprops(em_bw, im(:,:,1,t,ch), 'Centroid',...
%             'MajorAxisLength', 'MinorAxisLength', 'Orientation');
% 
%         if size(props_em,1) <= 1
%             major_axes = props_em.MajorAxisLength;
%             minor_axes = props_em.MinorAxisLength;
%     
%             % Save centers of nuclei
%             em_center = props_em.Centroid;
%             em_angle = props_em.Orientation;
%             xy = calc_ellipse_params(em_angle, major_axes, minor_axes,...
%                 em_center);
%             
%             ac = zeros(1,2);
%             ac(1) = (xy(2,2) - xy(1,2)) ./ (xy(2,1) - xy(1,1));
%             ac(2) = xy(2,2) - (ac(1) * xy(2,1));
%         end

        if (t ~= 1) && size(centers{t}, 1) < 200
            bw = mask(:,:,1,t-1);
        end

        mask(:,:,1,t) = bw;

        CC = bwconncomp(bw);
        L = labelmatrix(CC);

        props = regionprops(L, im(:,:,1,t,ch),...
                'Centroid', 'MeanIntensity', 'MajorAxisLength',...
                'MinorAxisLength', 'PixelIdxList', 'Circularity');
        
        % Save the mean intensities
        avg_I{t} = cat(1, props.MeanIntensity);

        % Save centers of nuclei
        centers{t} = cat(1, props.Centroid);

        ellipse_axes{t} = cat(2, cat(1, props.MajorAxisLength), ...
                             cat(1, props.MinorAxisLength));

        % Save list of pixels
        ind{t} = {props.PixelIdxList}';

        circularity{t} = cat(1, props.Circularity);
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

    im_LoG = edge(im,'log',0,4);
    im_LoG = imfill(im_LoG, 'holes');
    se = strel('disk', 5);
    bw = imopen(im_LoG, se);

    bw = bwareafilt(bw, [100, 1000]);

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
    end
end

function [em_bw, boundary] = segment_embryo(im,T)
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

    % Segment entire embryo using Otsu's method for thresholding
    % Gaussian blur with standard deviation specified
    im_blur = imgaussfilt(im, 10);
    % Segment entire embryo from blurred image using thresholding
    
    em_mask = imbinarize(im_blur, T); % T = 0.02
    se = strel('disk',3);
    em_bw = imclose(em_mask,se);
    em_bw = bwareaopen(em_bw, 10000);

    B = bwboundaries(em_bw, 'noholes');
    
%     % Get parameters of an ellipse, alpha is angle between major axis and x
%     % axis, a is half the length of the major axis, b is half the length of
%     % the minor axis and c is the center of the ellipse
%     stats = regionprops(em_bw, 'Centroid', 'MajorAxisLength',...
%                             'MinorAxisLength', 'Orientation');
%     
%     alpha = -stats.Orientation * pi./180;
%     a = stats.MajorAxisLength/2;
%     b = stats.MinorAxisLength/2;
%     c = stats.Centroid;

    if length(B) == 1
       boundary = B{1};
       boundary(boundary(:,2) == 1,:) = [];
       boundary(boundary(:,2) == (size(im,2)), :) = [];
       boundary(boundary(:,1) == 1,:) = [];
       boundary(boundary(:,1) == (size(im,2)), :) = [];
    else
        boundary = [];
    end
end

% function xy = calc_ellipse_params(angle_deg, major, minor,...
%     center)
% %CALC_ELLIPSE_PARAMS Determines values for an ellipse and plots.
% %   
% %   Input
% %       uiax: the handle to the axes in the gui
% %       angle_deg: angle of ellipse from x-axis in degrees
% %       major: major axis of ellipse
% %       minor: minor axis of ellipse
% %       center: center of ellipse
% %       len: length/width of square pixels in microns
% %       vis: the flag for determing visibility of ellipse
% % 
% %   Output
% %       h: the handle to the ellipse
% 
%     % Calculate the values for plotting an ellipse
%     alpha = -angle_deg .* pi./180;
%     a = major./2;
%     b = minor./2;
%     c = center;
%     phi = [0,pi];
%     
%     % Calculate points along ellipse and plot
%     xy = parameterized_ellipse(phi, a, b, alpha, c);
% end

function xy = parameterized_ellipse(phi, a, b, alpha, c)
%PARAMETERIZED_ELLIPSE Calculates points on an ellipse.
%   
%   Input
%       phi: parametric angle
%       a: major axis
%       b: minor axis
%       alpha: angle between the major axis and the x-axis
%       c: center point
% 
%   Output
%       xy: xy positions of points defined by parametric equations of an
%           ellipse
    
    % Initilize variable for saving x and y
    xy = zeros(size(phi,2), 2);
    
    % Parametric equation for an ellipse
    X = a .* cos(phi);
    Y = b .* sin(phi);
    xy(:,1) = X .* cos(alpha) - Y .* sin(alpha) + c(1);
    xy(:,2) = X .* sin(alpha) + Y .* cos(alpha) + c(2);
end

% function d = calc_dist_to_line(ac, xy)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
%     
%     a = ac(:,1);
%     b = 1;
%     c = ac(:,2);
%     x = xy(:,1);
%     y = xy(:,2);
% 
%     d = abs(a.*x + b.*y + c) ./ sqrt(a.^2 + b.^2);
% end

function d = calc_dist_to_line(p, centers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    x1 = p(1,1);
    y1 = p(1,2);
    x2 = p(2,1);
    y2 = p(2,2);
    x0 = centers(1,1);
    y0 = centers(1,2);

    d = abs((x2 - x1) .* (y1 - y0) - (x1 - x0) .* (y2 - y1)) ./...
                    sqrt((x2 - x1).^2 + (y2 - y1).^2);
    
%     if x2 ~= x1
%         m = (y2 - y1)/(x2 - x1);
%         b_0 = y2 - (m .* x2);
%         y3 = m .* x0 + b_0;
% 
%         if y0 < y3
%             d = -d;
%         end
%     elseif (x0 < x2) && (x0 < x1)
%         d = -d;
%     end
end

function af = fit_ellipse(a0, x, y, im_len)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    f = @(a) (((x-a(1)) * cos(a(5)) + (y - a(2)) * sin(a(5))).^2)/a(3).^2 ...
        + (((x-a(1)) * sin(a(5)) - (y-a(2)) * cos(a(5))).^2) / a(4).^2 - 1;

    options = optimset('Display','off');
    ub = [im_len, im_len, Inf, Inf, pi/2];
    lb = [1, 1, 0, 0, -pi/2];
    af = lsqnonlin(f, a0, lb, ub, options);
end       