function [data] = calc_gene_props(varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    % If first input is a structure
    if isstruct(varargin{1})
        % n is size of structure, first input is the image data, second is 
        % the data structure
        n = size(varargin{1}, 2);
        data = varargin{1};
    else
        % Use menu to select files
        [name, folder] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
                'Select the microscope images', 'Multiselect', 'on');
        % n is first input and make the data structures for output
        % Calculate n
        n = size(name, 2);
    
        data = struct('folder', cell(1,n),...
                    'name', [],...
                    'raw_image', [],...
                    'pixel_length', [],...
                    'gene_length', [],...
                    'circle_overlap', [],...
                    'end_pts', [],...
                    'phi', [],...
                    'phi_ind', [],...
                    'arc_length', [],...
                    'perimeter', [],...
                    'z_plane', [],...
                    'bckgrd_channel', [],...
                    'signal_channels', [],...
                    'threshold', [],...
                    'n_domains', []);
    end

    for i = 1:n
        if ~isstruct(varargin{1})
            % Open image/movie
            [data(i).folder, data(i).name, data(i).raw_image,...
                data(i).pixel_length] = open_img(name{i}, folder);
            
            data(i).z_plane = varargin{1};
            data(i).bckgrd_channel = varargin{2};
            data(i).signal_channels = varargin{3};
            data(i).threshold = varargin{4};
            data(i).n_domains = varargin{5};
        end

        z = data(i).z_plane;
        ch = data(i).bckgrd_channel;
        channels_to_process = data(i).signal_channels;

        im = data(i).raw_image(:,:,z,ch);

        im_LoG = edge(im, 'log', 0, 20);

%         figure, imshow(im_LoG);

        props = regionprops(im_LoG, im, 'FilledArea', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation','PixelList');
        [~,I] = maxk(cat(1,props.FilledArea), 2);
        props = props(I);
        

%         bw_ring = bwareaopen(im_LoG, 1000);
% 
% 
%         bw_solid = imfill(bw_ring, 'holes');
%         bw_solid = bwareaopen(bw_solid, 50000);
%         bw_ring(~bw_solid) = false;

%         props = regionprops(bw_ring, im, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation','PixelList');

        af = zeros(size(props, 1), 5);

        for j = 1:size(props,1)
            x = props(j).PixelList(:,1);
            y = props(j).PixelList(:,2);
            a0 = [props(j).Centroid(1), props(j).Centroid(2), props(j).MajorAxisLength,...
            props(j).MinorAxisLength, props(j).Orientation.* pi./180];
    
            af(j,:) = fit_ellipse(a0, x, y, size(im, 1));
        end

        [xy, data(i).phi, a, b, alpha] = calc_ellipse_params(af);
        xy_round = round(xy);
        data(i).perimeter = integral(@(x) integrand(x, a, b, alpha), -pi, pi);

        data(i).gene_length = cell(size(channels_to_process,2),1);
        data(i).arc_length = cell(size(channels_to_process,2),1);
        data(i).circle_overlap = cell(size(channels_to_process,2),1);
        data(i).end_pts = cell(size(channels_to_process,2),1);
        data(i).phi_ind = cell(size(channels_to_process,2,1));

        for c = 1:size(channels_to_process,2)
            im_blur = imgaussfilt(data(i).raw_image(:,:,z,channels_to_process(c)), 9);
            bw_thresh = imbinarize(im_blur, data(i).threshold(c));
            bw_domain = bwareaopen(bw_thresh, 100);
%             figure; imshow(bw_domain);
            
            props2 = regionprops(bw_domain, data(i).raw_image(:,:,z,channels_to_process(c)), 'FilledArea','PixelList');
            if ~isempty(props2)
                [~,I2] = maxk(cat(1,props2.FilledArea), data(i).n_domains(c));
                props2 = props2(I2);

%                 disp(data(i).gene_length(c))
%                 pause
                data(i).gene_length{c} = zeros(data(i).n_domains(c), 1);
                data(i).arc_length{c} = zeros(data(i).n_domains(c), 1);
                data(i).circle_overlap{c} = cell(size(props2,1), 1);
                data(i).end_pts{c} = cell(size(props2,1), 1);
                data(i).phi_ind{c} = cell(size(props2,1), 1);
    
                for m = 1:size(props2,1)
                    data(i).circle_overlap{c}{m} = false(size(xy_round,1),1);
    
                    for k = 1:size(xy_round,1)
                        data(i).circle_overlap{c}{m}(k) = any(sum(props2(m).PixelList == xy_round(k,:),2) == 2);
                    end

                    if any(data(i).circle_overlap{c}{m})
                        data(i).end_pts{c}{m} = cat(1,diff(data(i).circle_overlap{c}{m}),...
                            data(i).circle_overlap{c}{m}(1)-data(i).circle_overlap{c}{m}(end));
            
                        phi_ind1 = find(data(i).end_pts{c}{m}==1);
                        phi_ind2 = find(data(i).end_pts{c}{m}==-1);
            
                        if (size(phi_ind1,1) > 1)
                            if (phi_ind2(1) > phi_ind1(1))
                                phi_ind1 = phi_ind1(1);
                            else
                                phi_ind1 = phi_ind1(end);
                            end
                            phi_ind2 = phi_ind2(end);
                        end
            
                        data(i).phi_ind{c}{m} = cat(2,phi_ind1,phi_ind2);
        
        %                 if i == 8 && c == 3
        %                     disp(data(i).phi_ind{c}{n})
        %                     pause
        %                 end
                        
                        if data(i).phi(data(i).phi_ind{c}{m}(1,1)+1) > data(i).phi(data(i).phi_ind{c}{m}(1,2))
                            data(i).arc_length{c}(m) = integral(@(x) integrand(x, a, b, alpha),...
                                (data(i).phi( data(i).phi_ind{c}{m}(1,1)+1) - 2*pi), data(i).phi(data(i).phi_ind{c}{m}(1,2)));
                        else
                            data(i).arc_length{c}(m) = integral(@(x) integrand(x, a, b, alpha),...
                                data(i).phi(data(i).phi_ind{c}{m}(1,1)+1), data(i).phi(data(i).phi_ind{c}{m}(1,2)));
                        end
                        
                        data(i).gene_length{c}(m) = data(i).arc_length{c}(m)/data(i).perimeter;
                    else
                        data(i).arc_length{c}(m) = 0;
                    end
                end
            else
                data(i).arc_length{c} = 0;
            end
        end
    end
end

function [path, file_name, im, pix_len] = open_img(name, folder)
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
    file_ext = strsplit(name, '.');

    file_name = file_ext{1};

    % Use bioformats to read in file
    reader = bfGetReader(path);
    omeMeta = reader.getMetadataStore();

    % Save the size of X, Y, Z, T, and C
    X = omeMeta.getPixelsSizeX(0).getValue();
    Y = omeMeta.getPixelsSizeY(0).getValue();
    Z = omeMeta.getPixelsSizeZ(0).getValue();
%     T = omeMeta.getPixelsSizeT(0).getValue();
    C = omeMeta.getPixelsSizeC(0).getValue();
    
    % Allocate looped variable
    I = uint8(zeros(X,Y,C,Z));
    
    % For each z slice
    for z = 1:Z
        % For each channel
        for c = 1:C
            % Get the index and save the image
            i = reader.getIndex(z-1, c-1, 0)+1;
            I(:,:,c,z) = bfGetPlane(reader, i);
        end
    end
    
    % Close the open file
    reader.close()
    
    % Save the physical length of a pixel in µm
    xy_len = omeMeta.getPixelsPhysicalSizeX(0).value(...
                            ome.units.UNITS.MICROMETER);
    
    % Convert pixel length to a double
    xy_len = xy_len.doubleValue();

    pix_len = cat(2, xy_len, xy_len);

    % Reshape image data to match dimensions, X, Y, z, time, channels
    im = permute(I, [1,2,4,3]);
end


function af = fit_ellipse(a0, x, y, im_len)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    f = @(a) (((x-a(1)) * cos(a(5)) + (y - a(2)) * sin(a(5))).^2)/a(3).^2 ...
        + (((x-a(1)) * sin(a(5)) - (y-a(2)) * cos(a(5))).^2) / a(4).^2 - 1;

    options = optimset('Display','off');
    ub = [im_len, im_len, im_len, im_len, pi/2];
    lb = [1, 1, 0, 0, -pi/2];
    af = lsqnonlin(f, a0, lb, ub, options);
end 

function [xy,phi,a,b,alpha] = calc_ellipse_params(af)
%CALC_ELLIPSE_PARAMS Determines values for an ellipse and plots.
%   
%   Input
%       uiax: the handle to the axes in the gui
%       angle_deg: angle of ellipse from x-axis in degrees
%       major: major axis of ellipse
%       minor: minor axis of ellipse
%       center: center of ellipse
%       len: length/width of square pixels in microns
%       vis: the flag for determing visibility of ellipse
% 
%   Output
%       h: the handle to the ellipse

    af_mean = mean(af,1);

    % Calculate the values for plotting an ellipse
    alpha = af_mean(5);
    a = af_mean(3);
    b = af_mean(4);
    c = [af_mean(1),af_mean(2)];
    phi = 0:0.01:2*pi;
    
    % Calculate points along ellipse and plot
    xy = parameterized_ellipse(phi, a, b, alpha, c);
end

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

% A function to calculate arclength using elliptical integrals
function qwe = integrand(phi, a, b, alpha)
    dx = -a * sin(phi) * cos(alpha) - b * cos(phi) * sin(alpha);
    dy = -a * sin(phi) * sin(alpha) + b * cos(phi) * cos(alpha);
    
    qwe = sqrt(dx.^2 + dy.^2);
end