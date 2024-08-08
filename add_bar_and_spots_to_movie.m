function add_bar_and_spots_to_movie(im, data, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
%     [name,folder] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
%             'Select the microscope images', 'MultiSelect', 'off');
% 
%     path = fullfile(folder,name);
%     img = bfopen(path);
%     X = img{1,4}.getPixelsSizeX(0).getValue();
%     Y = img{1,4}.getPixelsSizeY(0).getValue();
%     Z = img{1,4}.getPixelsSizeZ(0).getValue();
%     T = img{1,4}.getPixelsSizeT(0).getValue();
%     C = img{1,4}.getPixelsSizeC(0).getValue();
% 
%     % Intialize looped variable
%     raw_t = zeros(size(img{1,1},1),1);
% 
%     % For each image plane
%     for i = 1:size(img{1,1},1)
%         % Try to get the time that elapsed during image aquisition
%         % If unable to this it means the last z stack is incomplete
%         try
%             raw_t(i,1) = img{1,4}.getPlaneDeltaT(0,...
%                                     i-1).value.doubleValue./60;
%             delete_last_t = false;
%         catch
%             delete_last_t = true;
%         end
%     end
%     
%     % Reshape time to match the dimensions of channel, z, and time
%     raw_t = reshape(raw_t,C,Z,T);
%     
%     % Times indicate when image aquisition finished, add preceding 0 to get
%     % start of each z-stack and thus each time point
%     t_lbl = [0;squeeze(raw_t(end,end,1:end-1))];
% 
%     avg_t = mean(diff(t_lbl));
% 
%     img = permute(reshape(cat(3, img{1,1}{:,1}), Y, X, C, Z, T),...
%                   [1,2,4,5,3]);
%     img = max(img, [], 3);
%     img = squeeze(img);
%     img = img(:,:,:,ch);
%     
%     % Delete last timepoint if z-stack is incomplete
%     if delete_last_t
%         t_lbl(end) = t_lbl(end-1) + avg_t;
%     end

%     J = imrotate(img, -90);
%     J = flip(J, 2);

    colors = {'#C00000'; '#ED7D31'; '#FFC000'; '#70AD47'};

    currentFolder = pwd;

    f = figure;
    ax = axes('Parent', f);

    for j = 1:size(im,2)
        file_name_parts = strsplit(im(j).name, '.');
        mkdir(currentFolder, file_name_parts{1});
        img = im(j).raw_image(:,:,:,:,im(j).signal_channels);

        if ~isnan(data(j).blue_light)
            if size(data(j).blue_light, 3) == 1
                t = data(j).blue_light(2,1):data(j).blue_light(2,2);
            elseif size(data(j).blue_light, 3) == 2
                t = [data(j).blue_light(2,1,1):data(j).blue_light(2,2,1),data(j).blue_light(2,1,2):data(j).blue_light(2,2,2)];
            end
        else
            t = -1;
        end
        t_lbl = data(j).time - data(j).time(data(j).t_align);
        
%         J = flip(J, 1);
%         x = 8:23;
%         y = 12:112;
    %     y = 1024-(800:1000);
    %     x = 5:20;
    %     y = 5:105;
        D = duration(minutes(t_lbl),'format','hh:mm:ss');
    %     zoom_y = 260:465;
    %     zoom_x = 110:620;
        
        if (nargin >= 6) && ~isempty(varargin{4})
            t_range = varargin{4}{j};
        else
            t_range = 1:size(img,4);
        end

        if (nargin >= 5) && ~isempty(varargin{3})
            c = colors{varargin{3}(j)};
        else
            c = '#C00000';
        end

        for i = t_range
            if (nargin >=9) && ~isempty(varargin{7}) && strcmp(varargin{7}, 'RGB')
                imshow(squeeze(im(j).RGB_image(:,:,i,:)), 'Parent', ax);
            else
                temp_im = imadjust(img(:,:,1,i), [0/65536,9000/65536]);% [0/65536,30000/65536] and [0/65536,6000/65536], [0/65536,3000/65536], [0/65536,9000/65536]
                imshow(temp_im, 'Parent', ax);
            end
            hold on;

            if isfield(data, 'rm_pts') && ~isempty(data(j).rm_pts) 
                if ~isempty(data(j).rm_pts{i,1})
                    plot(ax,data(j).rm_pts{i,1}(:,1), data(j).rm_pts{i,1}(:,2), ...
                    'o', 'MarkerSize', 10, 'MarkerEdgeColor', c, 'Linewidth', 1.5);

                    if (nargin >= 7) && ~isempty(varargin{5}) && varargin{5}(j)
                        % Call plot_ellipse for the updated data
                        if (i < data(j).nuc_cycle(1,2)) && (i > data(j).nuc_cycle(1,1))
                            index_t = 1;
                        elseif (i < data(j).nuc_cycle(2,2)) && (i > data(j).nuc_cycle(2,1))
                            index_t = 2;
                        elseif (i < data(j).nuc_cycle(3,2)) && (i > data(j).nuc_cycle(3,1))
                            index_t = 3;
                        else
                            index_t = [];
                        end
                        
                        if ~isempty(index_t)
                            k = data(j).hull_pts{1,index_t};
                            plot(ax,data(j).pts{1,index_t}(k,1), data(j).pts{1,index_t}(k,2), '--', 'color', 'white', 'Linewidth', 1);
                        end
                    end
                end
            elseif 0%isfield(data, 'centers') && ~isempty(data(j).centers{i,1})
                plot(ax,data(j).centers{i,1}(:,1), data(j).centers{i,1}(:,2), ...
                'o', 'MarkerSize', 10, 'MarkerEdgeColor', c, 'Linewidth', 1.5);
            else
                plot(ax,-1,-1, ...
                'o', 'MarkerSize', 10, 'MarkerEdgeColor', c, 'Linewidth', 1.5);
            end

            if (nargin >= 4) && ~isempty(varargin{2})
                if varargin{2}(j,1)
                    if varargin{2}(j,2)
%                         img = flip(img, 1);
                        set(ax, 'xdir', 'reverse');
                    end

%                     img = imrotate(img, 90);
                    camroll(ax,90);
                else
                    if varargin{2}(j,2)
%                     img = flip(img, 1);
                        set(ax, 'ydir', 'normal');
                    end
                end
            end

            if (nargin >= 3) && strcmp(varargin{1}, 'movie')
%                 if i > 20
%                     color = 'white';
%                 else
%                     color = 'black';
%                 end

                text(0.99, 0.035, sprintf('%s hrs:min:sec',char(D(i))),...
                    'FontSize', 12, 'Color', 'white', 'Units', 'normalized', 'HorizontalAlignment', 'right');

                if (nargin >= 8) && ~isempty(varargin{6})
                    text(0.01, 0.945, varargin{6}{j}, 'FontSize', 12,...
                        'Color', 'white', 'Interpreter', 'latex', 'Units', 'normalized');
%                     text(0.01, 0.965, varargin{6}{j}, 'FontSize', 12,...
%                         'Color', 'white', 'Interpreter', 'latex', 'Units', 'normalized');
                end

                if isfield(data, 'bleach_frames') && any(i == data(j).bleach_frames)
                    text(0.77, 0.9, 'BLEACH', 'FontSize', 12,...
                        'Color', 'white', 'Interpreter', 'latex', 'Units', 'normalized');
                end
                
                if any(i == t)
                    h = annotation(f,'rectangle', [.7,.906,.12,.02], 'FaceColor', [0,0.4470,0.7410],'EdgeColor', [0,0.4470,0.7410], 'LineWidth', 2);
                else
                    h = [];
                end
            else
                h = [];
            end

            test = getframe(ax);

%             J = imrotate(mask_im, 90);
%             J = flip(J, 1);
%             J = im2uint8(J);
            hold off;
            clf(ax,'reset');
            if ~isempty(h)
                delete(h);
            end
    %         mask_im = rot90(mask_im, 3);
    %         mask_im = flip(mask_im, 2);
%             J = cat(3, J(:,:,1), J(:,:,1), J(:,:,1));

             if (nargin >= 4) && ~isempty(varargin{2})
                if varargin{2}(j,1) && ~varargin{2}(j,2)
                    J = test.cdata(2:(end-1), 2:end, :);
                else
                    J = test.cdata(1:(end-1), 2:end, :);
                end
             end

%             if (nargin >= 3) && strcmp(varargin{1}, 'movie')
%                 
% %                 if any(i == t)
% %     %                 J(x,y,1) = 0;
% %     %                 J(x,y,2) = 65535 .* 0.4470;
% %     %                 J(x,y,3) = 65535 .* 0.7410;
% %                 
% %                     J(x,y,1) = 0;
% %                     J(x,y,2) = 255 .* 0.4470;
% %                     J(x,y,3) = 255 .* 0.7410;
% %                 end
% %                 
% %                 RGB = insertText(J, [300, 480], sprintf('%s hrs:min:sec',char(D(i))),...
% %                     'FontSize', 20, 'TextColor', 'white', 'BoxOpacity', 0);
% %                 
% %                 if (nargin >= 8) && ~isempty(varargin{6})
% %                     RGB = insertText(RGB, [250, 0], varargin{6}{j},...
% %                         'FontSize', 20, 'TextColor', 'white', 'BoxOpacity', 0);
% %                 end
% %             else
% %                 RGB = J;
%             end
            
            imwrite(J, sprintf('%s/%s/t%02d.tif',currentFolder,file_name_parts{1},i), 'compression', 'none');
        end
    end

    close(f);
end

% function add_bar_to_movie(t, ch, t_lbl)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
%     
% %     [name,folder] = uigetfile({'*.czi', 'CZI files (*.czi)'},...
% %             'Select the microscope images', 'MultiSelect', 'off');
% J = cell(2,1);
% for i = 1:2
%     [name,folder] = uigetfile({'*.avi', 'AVI files (*.avi)'},...
%             'Select the microscope images', 'MultiSelect', 'off');
% 
%     path = fullfile(folder,name);
%     img = bfopen(path);
%     X = img{1,4}.getPixelsSizeX(0).getValue();
%     Y = img{1,4}.getPixelsSizeY(0).getValue();
%     Z = img{1,4}.getPixelsSizeZ(0).getValue();
%     T = img{1,4}.getPixelsSizeT(0).getValue();
%     C = img{1,4}.getPixelsSizeC(0).getValue();
% 
% %     img = permute(reshape(cat(3, img{1,1}{:,1}), X, Y, C, Z, T),...
% %                   [1,2,4,5,3]);
% %     img = max(img, [], 3);
% %     img = squeeze(img);
% %     img = img(:,:,:,ch);
%     img = squeeze(permute(reshape(cat(3, img{1,1}{:,1}), Y, X, C, Z, T),...
%                   [1,2,3,4,5]));
% 
%     J{i,1} = imrotate(img,-125);
% end
% 
%     
% 
%     x = 15:45;
% %     y = 900:1000;
%     y = 1024-(800:1000);
%     k = 1;
%     zoom_rg = 300:1050;
%     D = duration(minutes(t_lbl),'format','mm:ss');
%     for i = 1:size(J{1,1},4)
%         mask_im = J{1,1}(zoom_rg,:,:,i);
% %         mask_im = imadjust(img(:,:,i), [53/65536,20634/65536]);
% %         mask_im = rot90(mask_im, 3);
% %         mask_im = cat(3, mask_im(:,:,1), mask_im(:,:,1), mask_im(:,:,1));
% 
%         if any(i == t)
% %             mask_im(x,y,1) = 0;
% %             mask_im(x,y,2) = 65535 .* 0.4470;
% %             mask_im(x,y,3) = 65535 .* 0.7410;
%         
%         mask_im(x,y,1) = 0;
%             mask_im(x,y,2) = 255 .* 0.4470;
%             mask_im(x,y,3) = 255 .* 0.7410;
%         end
%         
%         RGB = insertText(mask_im,[1130, 680], sprintf('%s min:sec',char(D(i))),...
%             'FontSize', 40, 'TextColor', 'white', 'BoxOpacity', 0);
%         
%         imwrite(RGB, sprintf('t%02d.tif',i), 'compression', 'none');
%         k = k + 1;
%     end
% 
%     for i = 1:3
%         empty_space = zeros(size(J{1,1},1), size(J{1,1},2), 3);
%         empty_space = empty_space(zoom_rg,:,:);
%         RGB = insertText(empty_space,[1130, 680], sprintf('%s min:sec',char(D(k))),...
%             'FontSize', 40, 'TextColor', 'white', 'BoxOpacity', 0);
%         imwrite(RGB, sprintf('t%02d.tif',k), 'compression', 'none');
%         k = k + 1;
%     end
%     
%     for i = 1:size(J{2,1},4)
%         mask_im = J{2,1}(zoom_rg,:,:,i);
% %         mask_im = imadjust(img(:,:,i), [53/65536,20634/65536]);
% %         mask_im = rot90(mask_im, 3);
% %         mask_im = cat(3, mask_im(:,:,1), mask_im(:,:,1), mask_im(:,:,1));
% 
%         if any(i == t)
% %             mask_im(x,y,1) = 0;
% %             mask_im(x,y,2) = 65535 .* 0.4470;
% %             mask_im(x,y,3) = 65535 .* 0.7410;
%         
%         mask_im(x,y,1) = 0;
%             mask_im(x,y,2) = 255 .* 0.4470;
%             mask_im(x,y,3) = 255 .* 0.7410;
%         end
%         RGB = insertText(mask_im,[1130, 680], sprintf('%s min:sec',char(D(k))),...
%             'FontSize', 40, 'TextColor', 'white', 'BoxOpacity', 0);
%         imwrite(RGB, sprintf('t%02d.tif',k), 'compression', 'none');
%         k = k + 1;
%     end
% end