function [img, data] = quantify_ms2(varargin)
%Quantify_MS2 Segments images for foci, specifically MS2/MCP spots.
% 
%   Inputs
%       varargin:
%           1): image structure returned from this function
%           2): data structure returned from this function
% 
%       OR
% 
%       varargin:
%           1): number of movies to be analyzed
%           2): cell array with condition names
%           3): cell array with indicies for plotting
%           4): cell array with array of thresholds
%           5): cell array with array of channel indicies. Each entry
%               matches an index for a channel with signal
%           6): cell array with array of background indices.
%           7): cell array with string for 3D or 2D
%       Example 1: (1, {'dark'}, {1:19}, {[0.01, 0.02, 0.03]},...
%                   {1}, {0}, {'2D'});
%       Example 2: (2, {'dark', 'light'}, {1:19, 2:20},...
%                   {[0.01, 0.02, 0.03], [0.01, 0.02, 0.03]},...
%                   {[2, 3], [2, 3]}, {1, 1}, {'3D', '3D'});
%       
%       Size of cell arrays must match number of movies
% 
%   Outputs
%       img: structure with following fields
%           folder: the folder path containing the opened file
%           name: part of the file name occuring before the first space
%           condition: name of experimental condition to group
%               similarly treated experiments
%           signal_channels: indices denoting which channels in image
%               contain signal
%           background_channel: index for channel in an image used for
%               segmenting a nucleus
%           raw_image: raw maximum z-projection of images
%           threshold: thresholds for segmenting MS2 foci
%           mask: mask of segmented MS2 foci
%           image_dims: '2D' or '3D' to denote if a z-projection should be
%               made or not
%       data: structure with following fields
%           folder: the folder path containing the opened file
%           name: part of the file name occuring before the first space
%           condition: name of experimental condition to group
%               similarly treated experiments
%           pixel_length: length of pixel in microns
%           time: the time for each z-projection, taken as the time of the
%               last z-plane in a z-stack
%           raw_time: the time to finish each z-plane
%           A: cell array of area (number of pixels * pixel_lenth^2) for
%               each spot, for each time point (row) and threshold (col)
%           avg_I: same as A except each entry contains the mean intensity
%               of each spot
%           max_I: same as A except each entry contains the max intensity
%               of each spot
%           sum_I: same as A except each entry contains the sum of all
%               intensities in each spot
%           centers: cell array containing the center location (pixels) for
%               each time point (row) and threshold (col)
%           spot_indicies: linear index of the image corresponding to
%               pixels in the segmented spots
%           ind: range of values to be plotted, example 1:19
%           nan_ind: indices for determining if data should be changed to
%               NaN, generally not used
%           blue_light: an array [z1,z2;t1,t2] that specifies the first
%               z-plane (z1,t1) where blue light was detected, and the last
%               z-plane (z2,t2) where blue light was detected.
% 
%   Overview
%       This function takes in either the number of movies to be analyzed,
%       or the data structures that are outputs of this function. If the
%       data structures are entered, no additional inputs are necessary. If
%       the number of movies to be analyzed are entered, then the
%       conditions, indicies, indicies for NaN, and thresholds should be
%       entered. Thresholds will likely need to be determined for each
%       image, and are on a scale of 0 to 1, see segment_ms2 below. If the
%       number of enteries in each cell array input do not match, an error
%       will be thrown.
    
    % Used in error message if size of inputs don't match
    error_msg = {'conditions'; 'indices'; 'indices for NaN';...
                 'thresholds'; 'signal_channels'; 'background_channel';...
                 'image_dims'};
    % Field names for inputs
    img_field = {'condition'; 'threshold'; 'signal_channels';...
                 'background_channel'; 'image_dims'};
    % Field names for inputs
    data_field = {'condition'; 'ind'};
%     % Field names for replacing data with NaN
%     field = {'A'; 'avg_I', 'max_I, sum_I'};
    
    % If first input is a structure
    if isstruct(varargin{1})
        % n is size of structure, first input is the image data, second is 
        % the data structure
        n = size(varargin{1}, 2);
        img = varargin{1};
        data = varargin{2};
    else
        % n is first input and make the data structures for output
        n = varargin{1};
        img = struct('folder', cell(1,n),...
                       'name', [],...
                       'condition', [],...
                       'signal_channels', 1,...
                       'background_channel', [],...
                       'raw_image', [],...
                       'threshold', [],...
                       'mask', [],...
                       'mask_nuc', [],...
                       'image_dims', []);
        data = struct('folder', cell(1,n),...
                       'name', [],...
                       'condition', [],...
                       'pixel_length', [],...
                       'time', [],...
                       'raw_time', [],...
                       'A', [],...
                       'avg_I', [],...
                       'max_I', [],...
                       'sum_I', [],...
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
                  
        % For each input
        for k = 2:(size(error_msg, 1) + 1)
            % If the number of arguments matches input number
            if nargin >= k
                % If the input size is greater than the number of movies
                if size(varargin{k},2) > n && ~isempty(varargin{k})
                    % Throw an error for that input
                    error('Too many %s', error_msg{k-1});
                % If the input size is less than the number of movies
                elseif size(varargin{k},2) < n && ~isempty(varargin{k})
                    % Throw an error for that input
                    error('Too few %s', error_msg{k-1});
                end
            end
        end
    end
    
    % Indices to read inputted data into the img and data structure
    ind1 = [2, 4, 5, 6, 7];
    ind2 = [2, 3];

    % For each movie/data entry
    for i = 1:n
        % If the first input is not a structure
        if ~isstruct(varargin{1})
            % For each input
            for j = 1:size(img_field, 1)
                % If the input is entered
                if nargin >= j && ~isempty(varargin{j})
                    % Read the input into the corresponding field
                    img(i).(img_field{j}) = varargin{ind1(j)}{i};
                end
            end

            for k = 1:size(data_field, 1)
                % If the input is entered
                if nargin >= k && ~isempty(varargin{k})
                    % Read the input into the corresponding field
                    data(i).(data_field{k}) = varargin{ind2(k)}{i};
                end
            end
            
            % Open image/movie
            [img(i).folder, img(i).name, img(i).raw_image,...
                data(i).pixel_length, data(i).time,...
                data(i).raw_time] = open_img(img(i).image_dims);

            % Copy folder and filename to data structure
            data(i).folder = img(i).folder;
            data(i).name = img(i).name;
        end
        
        % Segment MS2 spots
        [img(i).threshold, img(i).mask, img(i).mask_nuc,...
            data(i).A, data(i).avg_I, data(i).max_I,...
            data(i).sum_I, data(i).centers, data(i).spot_indices,...
            data(i).ind] = segment_ms2(img(i).raw_image,...
                                       img(i).threshold,...
                                       img(i).signal_channels,...
                                       img(i).background_channel,...
                                       img(i).image_dims,...
                                       data(i).ind,...
                                       data(i).t_align);
        
%         % If NaN indices are given
%         if ~isempty(data(i).nan_ind)
%             % For each field with data in it
%             for f = 1:size(field,1)
%                 % Set value of given indices to NaN
%                 ind = false(size(data(i).n_spots));
%                 ind(data(i).nan_ind,:,:) = true;
%                 data(i).(field{f,1})(ind) = NaN;
%             end
%         end
    end
    
    % If none of the images has a nuclear channel
    if all([img.background_channel] == 0)
        img = rmfield(img, 'mask_nuc');
    end
end

function [path, embryo_number, img, vox_len, t, raw_t] = open_img(dims)
%OPEN_IMG Open a czi with a z-stack, a time series, and channels
% 
%   Inputs
%       dims: '2D' or '3D' to determine if a z-projection is made or not
% 
%   Outputs
%       path: the folder path containing the opened file
%       embryo_number: part of the file name before the first space
%       img: raw images or maximum z-projection of images
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
%     img = bfopen(path);
    
%     % Save sizes of images in all dimensions, including time and color
%     % channels
%     X = img{1,4}.getPixelsSizeX(0).getValue();
%     Y = img{1,4}.getPixelsSizeY(0).getValue();
%     Z = img{1,4}.getPixelsSizeZ(0).getValue();
%     T = img{1,4}.getPixelsSizeT(0).getValue();
%     C = img{1,4}.getPixelsSizeC(0).getValue();
    
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
%     raw_t = zeros(size(img{1,1},1),1);
    
%     % For each image plane
%     for i = 1:size(img{1,1},1)
%         % Try to get the time that elapsed during image aquisition
%         % If unable to, then the last z stack is incomplete
%         try
%             % Get time bewteen each z slice
%             raw_t(i,1) = img{1,4}.getPlaneDeltaT(0,...
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
    img = permute(I2, [1,2,5,4,3]);
    
%     % Reshape image data to match dimensions, X, Y, z, time, channels
%     img = permute(reshape(cat(3, img{1,1}{:,1}), Y, X, C, Z, T),...
%                   [1,2,4,5,3]);
    % Code to open large images in smaller parts
%     img1 = permute(reshape(cat(3, img{1,1}{1:(250.*C*Z),1}), Y, X, C, Z, 250),...
%                   [1,2,4,5,3]);
%     img1 = max(img1, [], 3);
%     img2 = permute(reshape(cat(3, img{1,1}{((250.*C*Z)+1):end,1}), Y, X, C, Z, T-250),...
%                   [1,2,4,5,3]);
%     img2 = max(img2, [], 3);
% 
%     img = cat(4, img1,img2);

%     % If 2D is desired
%     if strcmp(dims, '2D')
%         % Make a maximum z-projection
%         img = max(img, [], 3);
%     end
    
    % Delete last timepoint if z-stack is incomplete
    if delete_last_t
        img = img(:,:,:,1:(end-1),:);
        t = t(1:(end-1),1);
    end
end

function [T, mask, mask_nuc, A, I,...
          max_I, sum_I, centers, spot_indices, ind] =...
                            segment_ms2(im, T, ch, bkgd, dims, ind, t_align)
%SEGMENT_MS2 Segments MS2/MCP foci.
% 
%   Input
%       im: the raw image or maximum z-projected image
%       T: the thresholds for segmenting MS2 spots. The thresholds should
%           be entered as a linear array, with as many thresholds included
%           as desired
%       ch: The channel index to determine which channels to look for spots
%           in
%       bkgd: The index of a background channel, specifcally expected to be
%           cell nuclei
%       dims: '2D' or '3D' to determine if images should be processed as 3D
%           or 2D
%       ind: the indicies for the range of time points for plotting and
%           determing the first time point for normalizing. If not entered,
%           the function defaults to using the full set of timepoints
% 
%   Output
%       T: thresholds used to segment the MS2 spots. Uses the defaults of
%           [0.04, 0.06, 0.08] if no thresholds were entered
%       mask: mask of the image segmentation for MS2 spots
%       mask_nuc: mask of the image segmentation for the background, which
%           are usually cell nuclei
%       n: number of MS2 spots detected
%       A: number of pixels for each MS2 spot detected
%       I: mean intensity for each MS2 spot detected
%       max_I: max intensity for each MS2 spot detected
%       sum_I: sum intensity for each MS2 spot detected
%       centers: cell array containing the center location (pixels) for
%           each time point (row) and threshold (col)
%       spot_indicies: linear index of the image corresponding to
%           pixels in the segmented spots
%       ind: the indicies for the range of time points for plotting and
%           determing the first time point for normalizing. If not entered,
%           the function defaults to using the full set of timepoints
% 
%   Overview
%       The function works for either 8 bit or 16 bit images. Thresholds
%       likely need to be chosen for each individual imaging/microscopy
%       setup. Thresholds are on a scale of 0 to 1. T = graythresh(J); can
%       be used to approximate a threshold by Otsu's method. Images are
%       processed by first subtracting the background using a gaussian
%       filter to blur the images, then images are segmented using a
%       threshold. Small objects are then removed (objects that are 1 pixel
%       or less). The areas and intensities for each object are saved.
%   
%   Notes
%       We found for our images GFP/RFP thresholds:
%       0.04/0.01 gets all spots, but detects false positives
%       0.06/0.02 gets the most spots without false positives
%       0.08/0.03 doesn't get false positives, but misses spots
    
    % If ind is NaN, meaning it wasn't an input
    if isnan(ind)
        % ind is the index of the full range of time points
        ind = 1:size(im,4);
    end
    
    % If data is 3D use inputs
    if strcmp(dims, '3D')
        % Use names from 3D regionprops
        field_dim = {'VoxelIdxList', 'Volume'};
    % Elseif data is 2D use inputs
    elseif strcmp(dims, '2D')
        % Use names from regionprops
        field_dim = {'PixelIdxList', 'Area'};
    end
    
    % Allocate looped variables
    mask = false(cat(2, size(im(:,:,:,:,1)), size(T,2), size(ch,2)));
    centers = cell(size(im,4), size(T,2), size(ch,2));
    spot_indices = cell(size(im,4), size(T,2), size(ch,2));
%     n = zeros(size(im,4), size(T,2), size(ch,2));
    A = cell(size(im,4), size(T,2), size(ch,2));
    I = cell(size(im,4), size(T,2), size(ch,2));
    max_I = cell(size(im,4), size(T,2), size(ch,2));
    sum_I = cell(size(im,4), size(T,2), size(ch,2));
    
    % If background index is greater than 0
    if bkgd > 0
        % Allocate array for storing mask of background/nuclear channel
        mask_nuc = false(size(im(:,:,:,:,1)));
    else
        % Else, leave mask empty
        mask_nuc = [];
    end
    
    % If thresholds are empty
    if isempty(T)
        % Use default thresholds
        T = [0.04, 0.06, 0.08]; %thresholds for MCP_GFP
%         T = [0.01, 0.02, 0.03]; %thresholds for MCP_RFP
    end

%     T_inc = (1+(0.4/size(im,4)):(0.4/size(im,4)):(1.4));

%     temp_im = zeros(size(im,1),size(im,2),size(im,3), 'uint16');
%     for q = 1:3
%         test_im = im(:,:,:,nc(q,1):nc(q,2),1);
%         temp_im(:,:,:,q) = uint16(max(test_im, [], 4));
%     end
    
    t_proj = max(im(:,:,:,:,ch(1)), [], 4);
    im_blur = imgaussfilt(t_proj, 100);

    if size(t_proj, 1) == 512
        em_mask = imbinarize(im_blur, 0.25);
    elseif size(t_proj, 1) == 1024
        em_mask = imbinarize(im_blur, 0.05);
    end

    se = strel('disk',5);
    em_bw = imclose(em_mask,se);
    em_bw = bwareaopen(em_bw, 100);

    T_inc = 0;

    % For each time point
    for t = 1:size(im,4)
        % If background index is greater than 0
        if bkgd > 0
            % Segment nuclei
            mask_nuc(:,:,:,t) = segment_nuclei(im, 0.02, t, bkgd); %0.01
        end
        
        % For each channel
        for c = 1:size(ch,2)
            % Perform blurring using a gaussian filter and then subtract
            % blurred image from original image. The size of the blurring
            % is dependent on the magnification of the objective used in
            % imaging.
            
            % If image is 3D
            if strcmp(dims, '3D')
                J = im(:,:,:,t,ch(c));

                % Background subtraction using a gaussian filter to blur
                im_bg_sub = J - imgaussfilt3(J, [10, 10, 1]);
    
                % Perform a gaussian blur
                B = imgaussfilt3(im_bg_sub, [5, 5, 0.5]);
            % Elseif image is 2D
            elseif strcmp(dims, '2D')

                special = false;

                if special == false
                    % Median filter to remove salt/pepper noise
                    J = medfilt2(im(:,:,:,t,ch(c)));

    %                 % Background subtraction using a gaussian filter to blur
    %                 im_bg_sub = J - imgaussfilt(J, 5);
                    
                    if size(J,1) == 512
                        % Background subtraction using a gaussian filter to blur
                        im_bg_sub = J - medfilt2(J, [5,5]);
                    elseif size(J,1) == 1024
                        % Background subtraction using a gaussian filter to blur
                        im_bg_sub = J - medfilt2(J, [10,10]);
                    end
                    
                    % Perform a gaussian blur of standard deviation 1
                    B = imgaussfilt(im_bg_sub, 2);
                elseif special == true
                    J = imgaussfilt(im(:,:,:,t,ch(c)), 1.5);
                    im_bg_sub = J - medfilt2(J, [25,25]);
                    B = im_bg_sub;
                end
            end
            
            % Code for processing each z slice individually
%             im_bg_sub = zeros(size(im,1), size(im,2), size(im,3));
%             B = zeros(size(im,1), size(im,2), size(im,3));
%             for z = 1:size(im,3)
%                 J = medfilt2(im(:,:,z,t,ch(c)));
%                 im_bg_sub(:,:,z) = J - imgaussfilt(J, 10);
%                 
%                 % Perform a gaussian blur of standard deviation 1
%                 B(:,:,z) = imgaussfilt(im_bg_sub(:,:,z), 5);
%             end

            % For each threshold
            for i = 1:size(T,2)
                if (t > (t_align + 40)) && (t_align ~= 1)
%                     T_inc = T_inc + 0.00025;
%                     0.000005 / T
%                     T_inc = 0.00025*log(t-(t_align + 40)); % sna
%                     T_inc = 0.001*log(t-(t_align + 40)); % sog
                      T_inc = 0*log(t-(t_align + 40)); % sog
                end
                
                % Make a mask by applying a threshold to the blurred image
                bw = imbinarize(B, T(c,i) + T_inc);
                
                % If a nuclear mask was made
                if ~isempty(mask_nuc)
                    % Any signal outside the nuclei is removed
                    bw(~mask_nuc(:,:,:,t)) = false;
                end

                bw(~em_bw) = false;

                if size(J,1) == 512
                    % Remove small objects that are 1 pixel
                    bw = bwareaopen(bw, 2);
                elseif size(J,1) == 1024
                    % Remove small objects that are 1 pixel
                    bw = bwareaopen(bw, 1); % used to be 5
                end

                bw = imclearborder(bw);

                masked_im = im(:,:,:,t,ch(c));
                masked_im(~bw) = 0;

                if size(J,1) == 512
                    bw_mask = imbinarize(masked_im, 0.3);
                elseif size(J,1) == 1024
                    bw_mask = imbinarize(masked_im, 0.075); % threshold was 0.2
                    se = strel('disk',5);
                    bw_mask = imclose(bw_mask,se);
                end

                CC = bwconncomp(bw);
                L = labelmatrix(CC);

                for k = 1:CC.NumObjects
                    if ~(any((L==k) & bw_mask, 'all'))
                        bw(L==k) = false;
                    end
                end

                
                % If image is 3D
                if strcmp(dims, '3D')
                    % Get properties of the segmented objects
                    props = table2struct(regionprops3(bw, im_bg_sub,...
                                  'Volume', 'MeanIntensity', 'Centroid',...
                                  'VoxelIdxList', 'MaxIntensity'));

                    % If fewer than three spots are detected and all the
                    % spots are smaller than 5 voxels
                    if size(props, 1) < 3 && all([props.Volume] < 5)
                        % Assume no real spots are detected by setting
                        % everything in mask to false
                        bw = false(size(bw));
                        
                        % Get the properties when there are no real spots
                        props = table2struct(regionprops3(bw, im_bg_sub,...
                                  'Volume', 'MeanIntensity', 'Centroid',...
                                  'VoxelIdxList', 'MaxIntensity'));
                    end
                % Elseif the image is 2D
                elseif strcmp(dims, '2D')
                    % Get properties of the segmented objects
                    props = regionprops(bw, im_bg_sub,...
                                    'Area', 'MeanIntensity', 'Centroid',...
                                    'PixelIdxList', 'MaxIntensity');

                    % If fewer than three spots are detected and all the
                    % spots are smaller than 5 pixels
                    if size(props, 1) < 3 && all([props.Area] < 5)
                        % Assume no real spots are detected by setting
                        % everything in mask to false
                        bw = false(size(bw));

                        % Get the properties when there are no real spots
                        props = regionprops(bw, im_bg_sub,...
                                    'Area', 'MeanIntensity', 'Centroid',...
                                    'PixelIdxList', 'MaxIntensity');
                    end
                end
                
                % Allocate space for sum of intensities
                sum_I{t,i,c} = zeros(size(props, 1), 1);
                
                % For each spot
                for j = 1:size(props, 1)
                    % Save the sum of intensities for each spot
                    sum_I{t,i,c}(j,1) = sum(...
                                    im_bg_sub(props(j,1).(field_dim{1})));
                end
                
                % Save the mask of the segmented spots
                mask(:,:,:,t,i,c) = bw;
                              
%                 % Save the number of objects/MS2 spots that were segmented              
%                 n(t,i,c) = size(props,1);
                
                % Save the number of pixels per MS2 spot detected
                A{t,i,c} = cat(1, props.(field_dim{2}));
                
                % Save the intensity per MS2 spot detected
                I{t,i,c} = cat(1, props.MeanIntensity);
                
                % Save the max intensity per MS2 spot detected
                max_I{t,i,c} = double(cat(1, props.MaxIntensity));
                
                % Save the centers of each spot
                centers{t,i,c} = cat(1, props.Centroid);
%                 if ~isempty(centers{t,i,c})
%                     centers{t,i,c} = centers{t,i,c}((centers{t,i,c}(:,2) > 5) & (centers{t,i,c}(:,2) < (size(bw,1)-4)), :);
%                 end
                
                % Save the list of indices for each spot
                spot_indices{t,i,c} = {props.(field_dim{1})}';
            end
        end
    end
end

function [mask] = segment_nuclei(img, T, t, ch)
%SEGMENT_MS2 Segments MS2/MCP foci.
% 
%   Input
%       img: the raw image or maximum z-projected image
%       T: the thresholds for segmenting MS2 spots. The thresholds should
%           be entered as a linear array, with as many thresholds included
%           as desired
%       t: the index for the time point
%       ch: the index for the color channel
% 
%   Output
%       mask: mask of the image segmentation for nuclei
% 
%   Overview
%       Images are processed by first performing background subtraction,
%       and then a gaussian filter is used to blur the images. Images are
%       segmented using a threshold. This should segment the background
%       channel, usually for cell nuclei
%   
%   Notes
%       We found for our images GFP/RFP thresholds:
%       0.04/0.01 gets all spots, but detects false positives
%       0.06/0.02 gets the most spots without false positives
%       0.08/0.03 doesn't get false positives, but misses spots
    
    % Intialize looped variables
%     props = cell(size(img,4), 1);
    
%     J = medfilt3(img(:,:,:,i,1));
    % Background subtraction
    im_bg_subtract = img(:,:,:,t,ch) - imgaussfilt3(img(:,:,:,t,ch),...
                                [50, 1, 0.5]);

    % Perform a gaussian blur
    B = imgaussfilt3(im_bg_subtract, [5, 5, 0.5]);
    
    % Segment using a threshold
    bw = imbinarize(B, T);
%     se = strel('cuboid', [5 5 1]);
%     bw = imopen(bw, se);

    % Remove small objects that are 1 pixels or less
    mask = bwareaopen(bw, 2);

%     % Get properties of the segmented objects
%     props{i} = table2struct(regionprops3(bw, im_bg_subtract,...
%                       'Volume', 'MeanIntensity', 'Centroid',...
%                       'VoxelIdxList'));
end