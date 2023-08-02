function display_mask_3D(data, varargin)
%DISPLAY_MASK Overlays mask and displays it
% 
%   Inputs
%       data: structure containing image and mask
%       varargin: structure containing points and hull information
%   
%   Outputs
%       None, will open a GUI to display images with overlayed masks.
% 
%   Overview
%       Takes a data structure and opens a GUI for displaying images and
%       creating overlays of the mask on the image, defined by field. The
%       function can display an overlay for each entry in the structure.
    
    % Set defaults
    [i, channels, colors, ch, color, c, z, t, T, B, n, em_vis,...
        im_vis, im_field, mask_field, im_type, im_ind, hull_vis,...
        pts_vis] = set_defaults(data);

    % Save image and mask
    [im, mask] = set_image_and_mask(data, im_field, mask_field, i, n);

    % Create a uifigure
    fig_pos = [0, 0, 1200, 750];
    f = uifigure('Name', 'Mask Viewer', 'NumberTitle', 'off',...
                 'Position', fig_pos, 'Visible', 'off');
    set(f, 'doublebuffer', 'off');
    movegui(f, 'center');
    set(f, 'Visible', 'on');
    
    % Create a grid layout with 5 rows and 6 columns. The first 4 rows
    % should be the same size, and the last row should be big enough to fit
    % the sliders
    g1 = uigridlayout(f,[5, 6]);
    g1.RowHeight = {'1x', '1x', '1x', '1x', 'fit'};
    g1.ColumnWidth = {'2.25x', '1x', '1x', '1x', '1x', '2x'};

    % Create a uiaxes to load the images into the uifigure
    uiax = uiaxes(g1, 'XLim', [0 size(im, 2)], 'YLim', [0 size(im, 1)],...
                  'XTick', [], 'YTick', [],...
                  'XTickLabel', [], 'YTickLabel', []);

    % Hide the axes so its an image display
    uiax.XAxis.Visible = false;
    uiax.YAxis.Visible = false;

    % Set the axis so its in row 1 through 4 and column 2 through 5 of the
    % grid layout
    uiax.Layout.Row = [1, 4];
    uiax.Layout.Column = [2, 5];
    
    % Create a nested gridlayout with 11 rows and 4 columns for the buttons
    g2 = uigridlayout(g1, [11, 4], 'Scrollable','on');

    % Set the new gidlayout in the second and fourth row and first column
    g2.Layout.Row = [2, 4];
    g2.Layout.Column = 1;

    % Set the first and last row, and first and last column to 1x so that
    % the buttons are centered. Set the rest of the rows to fit so the
    % label is readable. Set the width of the second column to a fixed 100
    % so that it doesn't change.
    g2.RowHeight = {'1x', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit',...
                    'fit', 'fit', 'fit', '1x'};
    g2.ColumnWidth = {'1x', 125, 100, '1x'};
    
    % Create a nested gridlayout with 3 rows and 5 columns for the sliders
    g3 = uigridlayout(g1,[3, 5]);

    % Set the new gidlayout in the fifth row and first through sixth column
    g3.Layout.Row = 5;
    g3.Layout.Column = [1 6];

    % Set all the row heights to fit so that the sliders fit, and set the
    % width of the label column as well as the width of the <> buttons
    % column, but allow the slider column to fill remaining column width
    g3.RowHeight = {'fit', 'fit', 'fit'};
    g3.ColumnWidth = {75, '1x', 40, 40, 40};

    % Initialize the figure with the first image
    [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask, ch,...
        c, B, T, data, i, em_vis, im_vis, n, color, hull_vis, pts_vis,...
        varargin{:});
    
    % If image is a z-stack
    if size(im,3) > 1
        % Create a text label for the z slider
        uilabel(g3, 'Text', 'z slice', 'HorizontalAlignment', 'right');
        
        % Create a slider for controlling which z plane to display, also
        % set the limits and tick marks
        sld_z = uislider(g3,...
                    'Limits', [1, size(im,3)],...
                    'MajorTicks', 1:size(im,3),...
                    'MinorTicks', [],...
                    'ValueChangingFcn',...
                                @(sld,event) update_z(sld, event),...
                    'ValueChangedFcn',...
                                @(sld,event) update_z_pos(sld, event));

        % Create a push button for decreasing z
        minus_z = uibutton(g3,'push',...
                       'Text', '<',...
                       'HorizontalAlignment', 'center',...
                       'VerticalAlignment', 'center',...
                       'FontSize', 20,...
                       'Enable', false,...
                       'ButtonPushedFcn',...
                           @(btn, events) minusZButtonPushed(btn, events));
    
        % Create a push button for increasing z
        plus_z = uibutton(g3,'push',...
                       'Text', '>',...
                       'HorizontalAlignment', 'center',...
                       'VerticalAlignment', 'center',...
                       'FontSize', 20,...
                       'ButtonPushedFcn',...
                            @(btn, events) plusZButtonPushed(btn, events));

        % Create a field for entering which z plane to display
        edt_z = uieditfield(g3, 'numeric',...
                          'HorizontalAlignment', 'center',...
                          'FontSize', 14,...
                          'Limits', [1, size(im,3)+0.5],...
                          'LowerLimitInclusive','on',...
                          'UpperLimitInclusive','on',...
                          'Value', 1,...
                          'ValueChangedFcn',...
                                    @(txt,event) edit_z_update(txt,event));
    end
    
    % If image is a time series
    if size(im,4) > 1
        tMinTick = ceil(size(im,4)/350) * 5;
        tMajTick = tMinTick * 5;

        % Create a text label for the time slider
        uilabel(g3, 'Text', 'time', 'HorizontalAlignment', 'right');
        
        % Create a slider for controlling which time point to display, also
        % set the limits and tick marks
        sld_t = uislider(g3,...
                    'Limits', [1, size(im,4)],...
                    'MajorTicks', [1, tMajTick:tMajTick:size(im,4)],...
                    'MinorTicks', [1, tMinTick:tMinTick:size(im,4)],...
                    'ValueChangingFcn',...
                                    @(sld,event) update_t(sld, event),...
                    'ValueChangedFcn',...
                                    @(sld,event) update_t_pos(sld, event));

        % Create a push button for decreasing time
        minus_t = uibutton(g3,'push',...
                       'Text', '<',...
                       'HorizontalAlignment', 'center',...
                       'VerticalAlignment', 'center',...
                       'FontSize', 20,...
                       'Enable', false,...
                       'ButtonPushedFcn',...
                           @(btn, events) minustButtonPushed(btn, events));
    
        % Create a push button for increasing time
        plus_t = uibutton(g3,'push',...
                       'Text', '>',...
                       'HorizontalAlignment', 'center',...
                       'VerticalAlignment', 'center',...
                       'FontSize', 20,...
                       'ButtonPushedFcn',...
                            @(btn, events) plustButtonPushed(btn, events));

        % Create a field for entering which time to display
        edt_t = uieditfield(g3, 'numeric',...
                          'HorizontalAlignment', 'center',...
                          'FontSize', 14,...
                          'Limits', [1 size(im,4)],...
                          'LowerLimitInclusive','on',...
                          'UpperLimitInclusive','on',...
                          'Value', 1,...
                          'ValueChangedFcn',...
                                    @(txt,event) edit_t_update(txt,event));
    end
    
    % If there is a threshold field and more than one threshold
    if isfield(data, 'threshold') && size(data(i).threshold,2) > 1
        % Create a text label for the threshold slider
        uilabel(g3, 'Text', 'Threshold', 'HorizontalAlignment', 'right');
        
        % Create a slider for controlling the displayed data for the
        % specified threshold and set the limits and tick marks
        sld_T = uislider(g3,...
                            'Limits', [1, size(data(i).threshold,2)],...
                            'MajorTicks', 1:size(mask,5),...
                            'MinorTicks', [],...
                            'ValueChangingFcn',...
                                    @(sld,event) update_T(sld, event),...
                            'ValueChangedFcn',...
                                    @(sld,event) update_T_pos(sld, event));

        % Create a push button for decreasing threshold
        minus_T = uibutton(g3,'push',...
                       'Text', '<',...
                       'HorizontalAlignment', 'center',...
                       'VerticalAlignment', 'center',...
                       'FontSize', 20,...
                       'Enable', false,...
                       'ButtonPushedFcn',...
                           @(btn, events) minusTButtonPushed(btn, events));
    
        % Create a push button for increasing threshold
        plus_T = uibutton(g3,'push',...
                       'Text', '>',...
                       'HorizontalAlignment', 'center',...
                       'VerticalAlignment', 'center',...
                       'FontSize', 20,...
                       'ButtonPushedFcn',...
                            @(btn, events) plusTButtonPushed(btn, events));

        % Create a entry field for specifying the threshold for the
        % displayed data
        edt_T = uieditfield(g3, 'numeric',...
                          'HorizontalAlignment', 'center',...
                          'FontSize', 14,...
                          'Limits', [1, size(data(i).threshold,2)],...
                          'LowerLimitInclusive', 'on',...
                          'UpperLimitInclusive', 'on',...
                          'Value', 1,...
                          'ValueChangedFcn',...
                                    @(txt,event) edit_T_update(txt,event));
    end
    
    % Create a text label for the brightness slider
    uilabel(g3, 'Text', 'Brightness', 'HorizontalAlignment', 'right');
    
    % Create a slider for controlling the brightness
    uislider(g3,...
                'Value', 100,...
                'Limits', [0, 100],...
                'MajorTicks', 0:5:100,...
                'ValueChangingFcn', @(sld,event) update_B(sld, event));
    
    % Label for the image selecting dropdown menu
    lbl_exp = uilabel(g2, 'Text', 'Select image',...
                          'HorizontalAlignment', 'left',...
                          'WordWrap', 'on');
    lbl_exp.Layout.Row = 2;
    lbl_exp.Layout.Column = 2;

    % Create a dropdown menu for changing which movie/image to display
    exp = uidropdown(g2, 'Items', {data.name},...
                         'Value', data(1).name,...
                         'ValueChangedFcn',...
                                      @(dd,event) exp_selection(dd,event));
    exp.Layout.Row = 2;
    exp.Layout.Column = 3;
    
    % Label for the mask overlaying button
    lbl_mask_btn = uilabel(g2, 'Text', 'Toggle on/off overlay of mask',...
                          'HorizontalAlignment', 'left',...
                          'WordWrap', 'on');
    lbl_mask_btn.Layout.Row = 3;
    lbl_mask_btn.Layout.Column = 2;

    % Create a push button for overlaying the mask on the image
    mask_btn = uibutton(g2,'state',...
                   'Text', 'Display Mask',...
                   'valueChangedFcn',...
                             @(mask_btn, events) maskButtonPushed(...
                                mask_btn, events));
    mask_btn.Layout.Row = 3;
    mask_btn.Layout.Column = 3;

    % Label for the circling detected spots button
    lbl_circle_btn = uilabel(g2, 'Text', 'Toggle on/off spots',...
                          'HorizontalAlignment', 'left',...
                          'WordWrap', 'on');
    lbl_circle_btn.Layout.Row = 4;
    lbl_circle_btn.Layout.Column = 2;

    % Create a push button for circling the detected spots on the image
    circle_btn = uibutton(g2,'state',...
                   'Text', 'Circle Spots',...
                   'valueChangedFcn',...
                             @(circle_btn, events) circleButtonPushed(...
                                circle_btn, events));
    circle_btn.Layout.Row = 4;
    circle_btn.Layout.Column = 3;

    % Label for the plotting convex hull button
    lbl_hull_btn = uilabel(g2, 'Text', 'Toggle on/off convex hull',...
                          'HorizontalAlignment', 'left',...
                          'WordWrap', 'on');
    lbl_hull_btn.Layout.Row = 5;
    lbl_hull_btn.Layout.Column = 2;

    % Create a push button for plotting the convex hull on the image
    hull_btn = uibutton(g2,'state',...
                   'Text', 'Display Convex Hull',...
                   'valueChangedFcn',...
                             @(hull_btn, events) hullButtonPushed(...
                                hull_btn, events));
    hull_btn.Layout.Row = 5;
    hull_btn.Layout.Column = 3;

    % Label for the dropdown menu for selecting which color channel to
    % display
    lbl_dd1 = uilabel(g2, 'Text', 'Choose image channel',...
                          'HorizontalAlignment', 'left',...
                          'WordWrap', 'on');
    lbl_dd1.Layout.Row = 6;
    lbl_dd1.Layout.Column = 2;
    
    % Create a dropdown menu for selecting which color channel to display
    dd1 = uidropdown(g2,...
                  'Items', channels,...
                  'Value', 'Channel 1',...
                  'ValueChangedFcn', @(dd,event) channel_select(dd,event));
    dd1.Layout.Row = 6;
    dd1.Layout.Column = 3;

    % Label for the dropdown menu to choose which color to display the mask
    % in
    lbl_dd2 = uilabel(g2, 'Text', 'Choose color of mask',...
                          'HorizontalAlignment', 'left',...
                          'WordWrap', 'on');
    lbl_dd2.Layout.Row = 7;
    lbl_dd2.Layout.Column = 2;

    % Create a dropdown menu to choose which color to display the mask in
    dd2 = uidropdown(g2,...
                  'Items', colors,...
                  'Value', 'Red',...
                  'ValueChangedFcn', @(dd,event) color_select(dd,event));
    dd2.Layout.Row = 7;
    dd2.Layout.Column = 3;

    % Label for the dropdown menu for selecting which mask to display
    lbl_dd3 = uilabel(g2, 'Text', 'Choose mask type',...
                          'HorizontalAlignment', 'left',...
                          'WordWrap', 'on');
    lbl_dd3.Layout.Row = 8;
    lbl_dd3.Layout.Column = 2;
    
    % Create a dropdown menu for selecting which mask to display
    dd3 = uidropdown(g2,...
                    'ValueChangedFcn', @(dd,event) mask_select(dd,event));
    dd3.Layout.Row = 8;
    dd3.Layout.Column = 3;
    
    % The following if statement determines if the data is from
    % quantify_ms2.m or from quantify_in_situ.m
    % If max_projections is a field
    if isfield(data, 'max_projection')
        % If size of field is 2
        if size(im_field, 1) == 2
            % Set the choices for the mask dropdown menu to differentiate
            % between a max projection and sum projections for both signal
            % and background
            dd3.Items = {'Max Signal'; 'Max Background';...
                        'Sum Signal'; 'Sum Background'};
            dd3.Value = 'Max Signal';
        % Else, set the choices for the mask dropdown menu to signal and
        % background
        else
            dd3.Items = {'Signal'; 'Background'};
            dd3.Value = 'Signal';
        end

        % Create a text label for the data set/image slider
        lbl_i = uilabel(g3, 'Text', 'Image',...
                            'HorizontalAlignment', 'right');
        lbl_i.Layout.Row = 2;
        lbl_i.Layout.Column = 1;

        % Create a slider for changing the data set/image
        sld_i = uislider(g3,...
                    'Limits', [1, size(data,2)],...
                    'MajorTicks', 1:ceil(size(data,2)/20):size(data,2),...
                    'MinorTicks', 1:size(data,2),...
                    'ValueChangingFcn',...
                                    @(sld,event) update_i(sld, event),...
                    'ValueChangedFcn',...
                                    @(sld,event) update_i_pos(sld, event));

        % Create a push button for decreasing the data set/image
        minus_i = uibutton(g3,'push',...
                       'Text', '<',...
                       'Enable', false,...
                       'ButtonPushedFcn',...
                           @(btn, events) minusiButtonPushed(btn, events));
    
        % Create a push button for increasing the data set/image
        plus_i = uibutton(g3,'push',...
                       'Text', '>',...
                       'ButtonPushedFcn',...
                            @(btn, events) plusiButtonPushed(btn, events));

        % Label for the button for displaying embryo ellipse
        lbl_b1 = uilabel(g2, 'Text',...
                                'Toggle on/off ellipse fit to embryo',...
                             'HorizontalAlignment', 'left',...
                             'WordWrap', 'on');
        lbl_b1.Layout.Row = 9;
        lbl_b1.Layout.Column = 2;

        % Create a push button for displaying a plot of an ellipse that is
        % calculated from the segmented embryo
        b1 = uibutton(g2,'state',...
                 'Text', 'Embryo Ellipse',...
                 'Position', [183, 70, 150, 22],...
                 'valueChangedFcn',...
                           @(btn, events) plotEmButtonPushed(btn, events));
        b1.Layout.Row = 9;
        b1.Layout.Column = 3;

        % Label for the button for displaying signal ellipse
        lbl_b2 = uilabel(g2, 'Text',...
                                'Toggle on/off ellipse fit to signal',...
                             'HorizontalAlignment', 'left',...
                             'WordWrap', 'on');
        lbl_b2.Layout.Row = 10;
        lbl_b2.Layout.Column = 2;

        % Create a push button for displaing a plot of an ellipse that is
        % calculated from the segmented signal
        b2 = uibutton(g2,'state',...
                 'Text', 'Signal Ellipse',...
                 'Position', [367, 70, 150, 22],...
                 'valueChangedFcn',...
                           @(btn, events) plotImButtonPushed(btn, events));
        b2.Layout.Row = 10;
        b2.Layout.Column = 3;
    else
        % If size of field.mask is equal to 2
        if size(mask_field, 1) == 2
            % Set mask dropdown menu choices to different MS2 or nuclear
            dd3.Items = {'Nuclear', 'MS21', 'MS22'};
            dd3.Value = 'Nuclear';
        % Elseif the size is only 1, set the choice to MS2
        elseif size(mask_field, 1) == 1
            dd3.Items = {'MS2'};
            dd3.Value = 'MS2';
        end
    end
    
    % Change the z slice value when the slider is moved
    function update_z(~, event)
        z = update_slider(event, size(im,3), false, plus_z, minus_z,...
                            sld_z, edt_z);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Change the z slice slider value so it is discrete
    function update_z_pos(sld, event)
        update_slider_pos(sld, event, false, size(im,3));
    end
    
    % Callback when a z is entered to change the z slice
    function edit_z_update(~,event)
        z = edit_x_update(event, size(im,3), false, plus_z, minus_z,...
                            sld_z, edt_z);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end

    % Create the function for the callback for when the minus z button is
    % pushed
    function minusZButtonPushed(~, ~)
        z = plusminusXButtonPushed(z, size(im,3), plus_z, minus_z,...
                                    sld_z, edt_z, 'minus');

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Create the function for the callback to the plus z button
    function plusZButtonPushed(~, ~)
        z = plusminusXButtonPushed(z, size(im,3), plus_z, minus_z,...
                                    sld_z, edt_z, 'plus');

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Change the time value when the slider is moved
    function update_t(~, event)
        t = update_slider(event, size(im,4), true, plus_t, minus_t,...
                            sld_t, edt_t);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Change the time slider value so it is discrete
    function update_t_pos(sld, event)
        update_slider_pos(sld, event, true, size(im,4));
    end
    
    % Callback when a t is entered to change the time
    function edit_t_update(~,event)
        t = edit_x_update(event, size(im,4), true, plus_t, minus_t,...
                            sld_t, edt_t);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end

    % Create the function for the callback to the minus button for time
    % index
    function minustButtonPushed(~, ~)
        t = plusminusXButtonPushed(t, size(im,4), plus_t, minus_t,...
                                    sld_t, edt_t, 'minus');

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Create the function for the callback to the plus button for time
    % index
    function plustButtonPushed(~, ~)
        t = plusminusXButtonPushed(t, size(im,4), plus_t, minus_t,...
                                    sld_t, edt_t, 'plus');

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end

    % Change the threshold value when the slider is moved
    function update_T(~, event)
        T = update_slider(event, size(data(i).threshold,2), false,...
                            plus_T, minus_T, sld_T, edt_T);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Change the threshold slider value so it is discrete
    function update_T_pos(sld, event)
        update_slider_pos(sld, event, false, size(data(i).threshold,2));
    end

    % Callback when a number is entered to change the threshold
    function edit_T_update(~,event)
        T = edit_x_update(event, size(data(i).threshold,2), false,...
            plus_T, minus_T, sld_T, edt_T);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end

    % Create the function for the callback to the minus button for
    % threshold
    function minusTButtonPushed(~, ~)
        T = plusminusXButtonPushed(T, size(data(i).threshold,2), plus_T,...
                                minus_T, sld_T, edt_T, 'minus');

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Create the function for the callback to the plus button for threshold
    function plusTButtonPushed(~, ~)
        T = plusminusXButtonPushed(T, size(data(i).threshold,2), plus_T,...
                                minus_T, sld_T, edt_T, 'plus');

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Function to update value of brightness for image
    function update_B(~, event)
        % Get new value for brightness factor from current position of
        % slider
        B = event.Value;

        % To prevent an error where max and min values are equal, add a
        % small number to brightness factor to make it nonzero if zero
        if B == 0
            B = B + 0.000000001;
        end

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Create callback for the dropdown menu to select the experiment
    function exp_selection(dd,~)
        % Find the index of the experiment with experiment name that
        % matches the value selected in the dropdown menu
        i = find(strcmp({data.name}, dd.Value));
        
        % If plus_i and minus_i exist, meaning that data is from
        % quantify_in_situ_.m
        if exist('plus_i', 'var') && exist('minus_i', 'var')
            % Check limits to determine if plus or minus button is active
            check_limits(i, size(data,2), plus_i, minus_i)
        end
        
        % Set the slider for the experiment slider to the new experiment
        % index
        sld_i.Value = i;

        % Set the image and mask to those determined by the new experiment
        % index
        [im, mask] = set_image_and_mask(data, im_field, mask_field, i, n);
        
        % If slide for time exists
        if exist('sld_t', 'var')
            % Change t back to 1 and reset slider
            t = change_exp(sld_t, size(im,4), plus_t, minus_t, edt_t,...
                           true);
        end
        
        % If slide for z exists
        if exist('sld_z', 'var')
            % Change z back to 1 and reset slider
            z = change_exp(sld_z, size(im,3), plus_z, minus_z, edt_z,...
                           false);
        end
        
        % If slide for threshold exists
        if exist('sld_T', 'var')
            % Change T back to 1 and reset slider
            T = change_exp(sld_T, size(data(i).threshold,2), plus_T,...
                           minus_T, edt_T, false);
        end

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Create the function for the callback when the display mask button is
    % pushed
    function maskButtonPushed(mask_btn,~)
        % If the button is not clicked
        if mask_btn.Value == false
            % Set the color to 0,0,0 so the mask is not seen
            c = [0, 0, 0];
        % Elseif the button is clicked
        elseif mask_btn.Value == true
            % Set the color to the color from the dropdown menu
            c = pick_color(color);
        end

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end

    % Create the function for the callback when the circle spots button is
    % pushed
    function circleButtonPushed(circle_btn,~)
        % If the button is not clicked
        if circle_btn.Value == false
            % Set the points visibility flag to false
            pts_vis = false;
        % Elseif the button is clicked
        elseif circle_btn.Value == true
            % Set the points visibility flag to true
            pts_vis = true;
        end

        % Set visibility handle to show or hide convex hull based on if the
        % button is pressed
        if ~isempty(pts_vis)
            set(p_h, 'Visible', pts_vis);
        end
    end

    % Create the function for the callback when the convex hull button is
    % pushed
    function hullButtonPushed(hull_btn,~)
        % If the plot convex hull function is not pressed
        if hull_btn.Value == false
            % Set the hull visibility flag to false
            hull_vis = false;
        % Elseif the plot convex function is pressed
        elseif hull_btn.Value == true
            % Set the hull visibility flag to true
            hull_vis = true;
        end
        
        % Set visibility handle to show or hide convex hull based on if the
        % button is pressed
        if ~isempty(hull_vis)
            set(hull_h, 'Visible', hull_vis);
        end
    end

    % Create callback when the channel is changed from the dropdown menu
    function channel_select(dd,~)
        % Save the channel number based on which channel is picked
        switch dd.Value
            case 'Channel 1'
                ch = 1;
            case 'Channel 2'
                ch = 2;
            case 'Channel 3'
                ch = 3;
            case 'Channel 4'
                ch = 4;
        end

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Create callback for when the color is channged from the dropdown menu
    function color_select(dd,~)
        % Set color to the selected value
        color = dd.Value;
        
        % If the display mask button is currently pressed
        if mask_btn.Value == true
            % Change the color of the mask to the new color picked from the
            % dropdown menu
            c = pick_color(color);
        end
        
        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end

    % Create callback for dropdown menu that picks the mask type
    function mask_select(dd,~)
        % Save values n, index into image and mask fields
        im_i = strcmp(im_type, dd.Value);
        n = im_ind(im_i, :);
%         switch dd.Value
%             case 'Nuclear'
%                 n = 1;
%                 m = 1;
%             case 'MS2'
%                 n = 2;
%                 m = 1;
%             case 'MS21'
%                 n = 2;
%                 m = 1;
%             case 'MS22'
%                 n = 2;
%                 m = 2;
%             case 'Signal'
%                 n = 1;
%                 m = 1;
%             case 'Background'
%                 n = 2;
%                 m = 1;
%             case 'Max Signal'
%                 n = 1;
%                 m = 1;
%             case 'Max Background'
%                 n = 2;
%                 m = 1;
%             case 'Sum Signal'
%                 n = 1;
%                 m = 1;
%             case 'Sum Background'
%                 n = 2;
%                 m = 1;
%         end
        
        % Set the image and mask to those determined by the selection
        [im, mask] = set_image_and_mask(data, im_field, mask_field, i, n);

        % Set the thresholds slider to the new thresholds and initialize it
        % to the first threshold
        if exist('sld_T', 'var')
            T = change_exp(sld_T, size(data(i).threshold,2), plus_T,...
                           minus_T, edt_T, false);
        end

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Change the experiment/image number value when the slider is moved
    function update_i(~, event)
        i = update_slider(event, size(data,2), false, plus_i, minus_i,...
                            sld_i, []);
    
        % Set the dropdown menu to the experiment that is currently
        % selected
        exp.Value = data(i).name;

        % Set the image and mask for the new experiment selected
        [im, mask] = set_image_and_mask(data, im_field, mask_field, i, n);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Change the experiment/image slider value so it is discrete
    function update_i_pos(sld, event)
        update_slider_pos(sld, event, false, size(data,2));
    end
    
    % Create the function for the callback to the minus button for
    % experiment/image number
    function minusiButtonPushed(~, ~)
        i = plusminusXButtonPushed(i, size(data,2), plus_i, minus_i,...
                                    sld_i, [], 'minus');

        % Set the dropdown menu to the experiment that is currently
        % selected
        exp.Value = data(i).name;
        
        % Set the image and mask to those determined by the new experiment
        % index
        [im, mask] = set_image_and_mask(data, im_field, mask_field, i, n);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end
    
    % Create the function for the callback to the plus button for
    % experiment/image number
    function plusiButtonPushed(~, ~)
        i = plusminusXButtonPushed(i, size(data,2), plus_i, minus_i,...
                                    sld_i, [], 'plus');

        % Set the dropdown menu for choosing the experiment to the new
        % experiment
        exp.Value = data(i).name;

        % Set the image and mask to those determined by the new experiment
        % index
        [im, mask] = set_image_and_mask(data, im_field, mask_field, i, n);

        % Make the mask overlay and update the image
        [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin{:});
    end

    % Create the function for the callback to the plot embryo ellipse
    % function
    function plotEmButtonPushed(btn,~)
        % If the plot embryo ellipse function is not pressed
        if btn.Value == false
            % Set the ellipse visibility flag to false
            em_vis = false;
        % Elseif the plot embryo ellipse function is pressed
        elseif btn.Value == true
            % Set the ellipse visibility flag to true
            em_vis = true;
        end
        
        % Set visibility handle to show or hide ellipse based on if button
        % is pressed
        if ~isempty(em_h)
            set(em_h, 'Visible', em_vis);
        end
    end

    % Create the function for the callback to the plot image ellipse
    % function
    function plotImButtonPushed(btn,~)
        % If the plot embryo ellipse function is not pressed
        if btn.Value == false
            % Set the ellipse visibility flag to false
            im_vis = false;
        % Elseif the plot embryo ellipse function is pressed
        elseif btn.Value == true
            % Set the ellipse visibility flag to true
            im_vis = true;
        end

        % Set visibility handle to show or hide ellipse based on if button
        % is pressed
        if ~isempty(im_h)
            set(im_h, 'Visible', im_vis);
        end
    end
end

function [em_h, im_h, hull_h, p_h] = display_image(uiax, z, t, im, mask,...
                ch, c, B, T, data, i, em_vis, im_vis, n, color,...
                hull_vis, pts_vis, varargin)
%DISPLAY_IMAGE Overlays mask or plots points and convex hull and displays
% 
%   Inputs
%       uiax: axis handle to axis in the uifigure
%       z: z slice index
%       t: time index
%       im: image data
%       mask: mask data
%       ch: channel index
%       c: color array for overlaying mask
%       B: brightness factor
%       T: threshold index
%       data: data structure
%       i: experiment/image index
%       em_vis: true or false to determine if ellipse for the background/
%           embryo should be displayed
%       im_vis:true or false to determine if ellipse for the signal/
%           image should be displayed
%       n: the index to the image field, mask field, and mask dimension
%       color: the current color for the points
%       hull_vis: true or false to determine if convex hull is visible
%       pts_vis: true or false to determine if points are visible
%       varargin: structure containing points and hull information
%   
%   Outputs
%       em_h: handle to the background/embryo ellipse plot
%       im_h: handle to the signal/image ellipse plot
%       hull_h: handle to the convex hull plot
%       p_h: handle to the points plot
% 
%   Overview
%       The image is brightend using the set brightness. If a mask is
%       specified, it is overlaid on the image specified by the inputs. The
%       grayscale image is converted to RGB, and then the color is
%       determined by decreasing the appropriate channels to zero where the
%       mask is. Also plots the centroids and convex hull if the
%       corresponding buttons were pressed.
    
    % Clear the axes
    cla(uiax);
    
    % Set limits to size of image, so only image is displayed
    uiax.XLim = [0 size(im, 2)];
    uiax.YLim = [0 size(im, 1)];

    % Set the image based on the brightness factor B
    brighter = imadjust(im(:,:,z,t,ch), [0, B ./ 100]);
    
    % If the mask variable is not empty
    if ~isempty(mask)
        % If the image is 8 bit
        if isa(im, 'uint8') && ~isempty(mask)
            % Make a RGB image where the original image is grayscale
            mask_im = cat(3, uint8(~mask(:,:,z,t,T)).^c(1) .* brighter,...
                             uint8(~mask(:,:,z,t,T)).^c(2) .* brighter,...
                             uint8(~mask(:,:,z,t,T)).^c(3) .* brighter);
        % Elseif the image is 16 bit
        elseif isa(im, 'uint16')
            % Make a RGB image where the original image is grayscale
            mask_im = cat(3, uint16(~mask(:,:,z,t,T)).^c(1) .* brighter,...
                             uint16(~mask(:,:,z,t,T)).^c(2) .* brighter,...
                             uint16(~mask(:,:,z,t,T)).^c(3) .* brighter);
        else
            im = im(:,:,z,t,ch);
            % Normalize image by maximum so that it displays. The 0.5
            % factor is to approximate the brightness of original images
            im = im./max(im(:)) .* 0.5;
            brighter = imadjust(im, [0, B ./ 100]);
            mask_im = cat(3, double(~mask(:,:,z,t,T)).^c(1) .* brighter,...
                             double(~mask(:,:,z,t,T)).^c(2) .* brighter,...
                             double(~mask(:,:,z,t,T)).^c(3) .* brighter);
        end
    
        % If the data includes the field angle_em_max, which means the data
        % is from quanitfy_in_situ.m and can plot an ellipse
        if isfield(data, 'angle_em_max')
            % Call plot_ellipse for the updated data
            [em_h, im_h] = plot_ellipse(uiax, data(i), em_vis, im_vis, n);
        else
            em_h = [];
            im_h = [];
        end

        % If postion data is provided, determined by an extra input
        % argument and varargin is not empty
        if nargin >= 18 && ~isempty(varargin{1})
            % If rm_pts is a field and is not empty
            if isfield(varargin{1}, 'rm_pts')...
                    && ~isempty(varargin{1}(i).rm_pts)...
                    && ~isempty(varargin{1}(i).rm_pts{t,1})
                hold(uiax,'on');
                % Plot points from rm_pts which is the list of points with
                % outliers removed
                p_h = plot(uiax, varargin{1}(i).rm_pts{t,1}(:,1),...
                    varargin{1}(i).rm_pts{t,1}(:,2), 'o',...
                    'MarkerSize', 12,...
                    'MarkerEdgeColor', color,...
                    'Linewidth', 1,...
                    'Visible', pts_vis);
                hold(uiax,'off');
            % Elseif centers is a field and is not empty
            elseif isfield(varargin{1}, 'centers')...
                    && ~isempty(varargin{1}(i).centers{t,1})
                hold(uiax,'on');
                % Plot points from centers which is the list of points
                % without outliers removed
                p_h = plot(uiax, varargin{1}(i).centers{t,1}(:,1),...
                    varargin{1}(i).centers{t,1}(:,2), 'o',...
                    'MarkerSize', 12,...
                    'MarkerEdgeColor', color,...
                    'Linewidth', 1,...
                    'Visible', pts_vis);
                hold(uiax,'off');
            % Else
            else
                % set p_h to empty
                p_h = [];
            end

            % If the data includes the field nuc_cycle, which means the
            % data is capable of plotting convex hull, and the pts and
            % hull_pts are not empty
            if isfield(varargin{1}, 'nuc_cycle')...
                    && ~isempty(varargin{1}(i).pts)...
                    && ~isempty(varargin{1}(i).hull_pts)
                % If t is less than the bounds of the nuclear cycle
                % timepoints, determine which index t is a part of and set
                % true false flag to true for the corresponding nuclear
                % cycle
                if (t < varargin{1}(i).nuc_cycle(1,2))...
                        && (t > varargin{1}(i).nuc_cycle(1,1))
                    index_t = 1;
                    tf_t = true;
                elseif (t < varargin{1}(i).nuc_cycle(2,2))...
                        && (t > varargin{1}(i).nuc_cycle(2,1))
                    index_t = 2;
                    tf_t = true;
                elseif (t < varargin{1}(i).nuc_cycle(3,2))...
                        && (t > varargin{1}(i).nuc_cycle(3,1))
                    index_t = 3;
                    tf_t = true;
                % Else
                else
                    % Set true false flag to false
                    tf_t = false;
                end
                
                % If true false flag is true
                if tf_t
                    % Save the indicis of the points that form the hull
                    k = varargin{1}(i).hull_pts{1,index_t};
                    hold(uiax,'on');

                    % Plot the points of the convex hull with lines
                    % connecting them
                    hull_h = plot(uiax,...
                        varargin{1}(i).pts{1,index_t}(k,1),...
                        varargin{1}(i).pts{1,index_t}(k,2), 'white',...
                        'Visible', hull_vis);
                    hold(uiax,'off');
                % Else if not plotting hull
                else
                    hull_h = [];
                end
            % Else if not plotting hull
            else
                hull_h = [];
            end
        % Else if not plotting points or hull
        else
            p_h = [];
            hull_h = [];
        end
    % Else if not diplaying mask, points, hull, or ellipses
    else
        mask_im = brighter;
        em_h = [];
        im_h = [];
        p_h = [];
        hull_h = [];
    end
    
    % Convert to 8 bit and display
    mask_im = im2uint8(mask_im);
    hold(uiax,'on');
    imshow(mask_im, 'Parent', uiax);
    hold(uiax,'off');
    
    % Set order of the handles of plotted objects to visualize on the image
    uistack(em_h,'top');
    uistack(im_h,'top');
    uistack(hull_h,'top');
    uistack(p_h,'top');
end

function [i, channels, colors, ch, color, c, z, t, T, B, n, em_vis,...
    im_vis, im_field, mask_field, im_type, im_ind, hull_vis,...
    pts_vis] = set_defaults(data)
%SET_DEFAULTS Sets values
% 
%   Inputs
%       data: data structure containing the image data
%   
%   Outputs
%       i: the experiment index
%       channels: the names of the channels
%       colors: the available colors for the mask
%       ch: the channel index
%       color: the current color for the mask
%       c: the color array, initially set so no mask is visible
%       z: the z stack index
%       t: the time index
%       T: the threshold index
%       B: the brightness factor, initially set so no brightening occurs
%       n: the index to the image field, mask field, and mask dimension
%       em_vis: the flag for if the plotted ellipse for the
%           background/embyo is visible
%       im_vis: the flag for if the plotted ellipse for the signal/image is
%           visible
%       im_field: a structure that contains the fieldnames for determining
%           which fields contains the imaging data
%       mask_field: a structure that contains the fieldnames for
%           determining which fields contain the mask data
%       im_type: an array with the names of possible names for the drop
%           down menu
%       im_ind: an array with the indices that correspond with the
%           image/mask selection
%       hull_vis: true/false to determine if convex hull plot is visible
%       pts_vis: true/false to determine if points plot is visible
% 
%   Overview
%       This function sets the parameters for the GUI.
    
    i = 1;
    colors = {'Red','Green','Blue','Cyan','Magenta','Yellow'};
    ch = 1;
    color = 'Red';
    c = [0, 0, 0];
    z = 1;
    t = 1;
    T = 1;
    B = 100;
    n = [1,1,1];
    em_vis = false;
    im_vis = false;
    hull_vis = false;
    pts_vis = false;
    
    % Save all field names
    field_names = fieldnames(data);

    % Find the field names that match image or projection but doesn't
    % include _dims. cellfun converts isempty to work on cell arrays to
    % give a true/false where the field name matched
    im_field = field_names(~cellfun('isempty',...
                                    regexp(field_names,...
                                            'image(?!_dims)|projection')));

    % Find the field names that match mask. cellfun converts isempty to
    % work on cell arrays to give a true/false where the field name matched
    mask_field = field_names(~cellfun('isempty',...
                                      regexp(field_names, 'mask')));
    
    % An array with the names of possible names for the drop down menu
    im_type = {'Nuclear'; 'MS2'; 'MS21'; 'MS22'; 'Signal'; 'Background';...
           'Max Signal'; 'Max Background'; 'Sum Signal'; 'Sum Background'};
    
    % An array with the indices that correspond with the image/mask
    % selection
    im_ind = [1,1,1; 1,2,1; 1,2,1; 1,2,2; 1,1,1; 1,2,1; 1,1,1; 1,2,1;...
              2,3,1; 2,4,1];

    % Save image using first field name
    im = data(i).(im_field{1});
    
%     % If the data is produced from quantify_in_situ rearrange the
%     % dimensions
%     if isfield(data, 'max_projection')
%         im = permute(im, [1,2,4,5,3]);
%     end

    % Create an array for channel names based on number of channels
    channels = append('Channel ', string(1:size(im, 5)));
end

function [im, mask] = set_image_and_mask(data, im_field, mask_field, i, n)
%SET_IMAGE_AND_MASK Updates the image and mask
% 
%   Inputs
%       data: data structure with images
%       im_field: structure with fieldnames corresponding to the correct
%           fields in data
%       i: experiment/image index
%       n: the index to the image field, mask field, and mask dimension
%   
%   Outputs
%       im: image corresponding to selected settings
%       mask: mask corresponding to selected settings
% 
%   Overview
%       Determines which image and which mask to display based on user
%       selection
    
    % Save image using field name that matches the index in n(1)
    im = data(i).(im_field{n(1)});
    
%     % If the data is produced from quantify_in_situ rearrange the
%     % dimensions
%     if isfield(data, 'max_projection')
%         im = permute(im, [1,2,4,5,3]);
%     end
    
    % Save mask using mask field name that matches index in n(2) or for the
    % mask in n(3)
    mask = data(i).(mask_field{n(2)})(:,:,:,:,:,n(3));
end

function c = pick_color(color)
%PICK_COLOR Picks color for the mask
% 
%   Inputs
%       color: color for mask, inserted as a string
%   
%   Outputs
%       c: array corresponding to color
% 
%   Overview
%       Takes a string (example: 'red') and uses it to pick the color of
%       the mask for displaying an overlay mask. The array c corresponds to
%       the exponents for determing the color of the mask by setting the
%       mask to zero in those channels/dim. Thus array c is the inverse of
%       a typical color array.
    
    % Corresponding array for matching color
    switch color
        case 'Red'
            c = [0, 1, 1];
        case 'Green'
            c = [1, 0, 1];
        case 'Blue'
            c = [1, 1, 0];
        case 'Cyan'
            c = [1, 0, 0];
        case 'Magenta'
            c = [0, 1, 0];
        case 'Yellow'
            c = [0, 0, 1];
    end
end

function [em_h, im_h] = plot_ellipse(uiax, data, em_vis, im_vis, n)
%PLOT_ELLIPSE Determines which values to use for plotting an ellipse.
%   
%   Input
%       uiax: the handle to the axes in the gui
%       data: the structure containing the data
%       em_vis: the flag for determing visibility of background/embryo
%           ellipse
%       im_vis: the flag for determing visibility of signal/image ellipse
%       n: the index for determing type of image/projection
% 
%   Output
%       em_h: the handle to the background/embryo ellipse
%       im_h: the handle to the signal/image ellipse
    
    % If using max projection
    if n(1) == 1
        % Calculate the values for plotting an ellipse using embryo or
        % background data for max projection
        em_h = calc_ellipse_params(uiax, data.angle_em_max,...
            data.Maj_em_max, data.Min_em_max, data.C_em_max,...
            data.pix_len, em_vis, 'w');
    
        % Calculate the values for plotting an ellipse using signal for max
        % projection
        im_h = calc_ellipse_params(uiax, data.angle_im_max,...
            data.Maj_im_max, data.Min_im_max, data.C_im_max,...
            data.pix_len, im_vis, 'y');
    % Else if using sum projection
    elseif n(1) == 2
        % Calculate the values for plotting an ellipse using embryo or
        % background data for sum projection
        em_h = calc_ellipse_params(uiax, data.angle_em_sum,...
            data.Maj_em_sum, data.Min_em_sum, data.C_em_sum,...
            data.pix_len, em_vis, 'w');
        
        % Calculate the values for plotting an ellipse using signal for sum
        % projection
        im_h = calc_ellipse_params(uiax, data.angle_im_sum,...
            data.Maj_im_sum, data.Min_im_sum, data.C_im_sum,...
            data.pix_len, im_vis, 'y');
    end
end

function h = calc_ellipse_params(uiax, angle_deg, major, minor,...
    center, len, vis, color)
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

    % Calculate the values for plotting an ellipse
    alpha = -angle_deg .* pi./180;
    a = major./(2 .* len);
    b = minor./(2 .* len);
    c = center;
    phi = 0:0.01:2*pi;
    
    % Calculate points along ellipse and plot
    xy = parameterized_ellipse(phi, a, b, alpha, c);
    hold(uiax,'on');
    h = plot(uiax, xy(:,1), xy(:,2), color, 'LineWidth', 3,...
        'Visible', vis);
    hold(uiax,'off');
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

function x = update_slider(event, upper_lim, val_flag, plus_x, minus_x,...
                            sld_x, edt_x)
%UPDATE_SLIDER Updates a value based on a changing slider
%   
%   Input
%       event: event where slider x moves
%       upper_lim: maximum or upper limit for slider value
%       val_flag: true or false for determining how to calculate slider
%           values
%       plus_x: uibutton handle for plus button for x
%       minus_x: uibutton handle for minus button for x
%       sld_x: slider handle for x
%       edt_x: editable field handle for x
% 
%   Output
%       x: updated value for x
%
%   Overview
%       Updates the value of x when the corresponding x slider is changing.
%       In addition, the function makes sure the corresponding plus and
%       minus buttons are available depending on if the max or min is
%       reached. The value in an editable field is also updated. This
%       function works generically for any slider (z, t, T etc).
    
    % If val_flag is true
    if val_flag
        % Make a temp array of all values
        tick_temp = 1:upper_lim;

        % Find the discrete value closest to the slider position
        [~, minIdx] = min(abs(event.Value - tick_temp(:)));

        % Override the selected value to the closest value
        x = tick_temp(minIdx);
    else
        % determine which discrete option the current value is closest
        % to for the slider when the slider is being moved, using the major
        % tick
        [~, minIdx] = min(abs(event.Value - event.Source.MajorTicks(:)));

        % Override the selected value
        x = event.Source.MajorTicks(minIdx);
    end
    
    % Check the current position against the limits to determine if plus
    % and minus buttons are clickable
    check_limits(x, upper_lim, plus_x, minus_x)
    
    % Set the value of the slider to the current value
    sld_x.Value = x;
    
    % If an editable field exists
    if ~isempty(edt_x)
        % Set the value in the editable field to the currently selected
        % value
        edt_x.Value = x;
    end
end

function update_slider_pos(sld, event, val_flag, upper_lim)
%UPDATE_SLIDER_POS Updates the slider position so it is discrete
%   
%   Input
%       sld: the slider handle
%       event: event where slider x moves
%       val_flag: true or false for determining how to calculate slider
%           values
%       upper_lim: maximum or upper limit for slider value
% 
%   Output
%       None.
%
%   Overview
%       This function moves the slider to discrete values so it is clear
%       that the slider cannot be intermmediate values. This function works
%       generically for any slider (z, t, T etc).
    
    % If val_flag is true
    if val_flag
        % Make a temp array of all values
        tick_temp = 1:upper_lim;

        % Find the discrete value closest to the slider position
        [~, minIdx] = min(abs(sld.Value - tick_temp(:)));

        % Move the slider to that option
        event.Source.Value = tick_temp(minIdx);
    else
        % determine which discrete option the current value is closest
        % to for the slider when the slider is moved, using the major tick
        [~, minIdx] = min(abs(sld.Value - event.Source.MajorTicks(:)));

        % Move the slider to that option
        event.Source.Value = event.Source.MajorTicks(minIdx);
    end
end

function x = edit_x_update(event, upper_lim, val_flag, plus_x, minus_x,...
                            sld_x, edt_x)
%EDIT_X_UPDATE Updates a value based on entered value
%   
%   Input
%       event: event where value for x is entered
%       upper_lim: maximum or upper limit for x value
%       val_flag: true or false for determining how to calculate values
%       plus_x: uibutton handle for > button for x
%       minus_x: uibutton handle for < button for x
%       sld_x: slider handle for x
%       edt_x: editable field handle for x
% 
%   Output
%       x: updated value for x
%
%   Overview
%       Updates the value of x when a value is entered into an editable
%       field for x. In addition, the function makes sure the corresponding
%       > and < buttons are available depending on if the max or min
%       is reached. The slider is also updated. This function works
%       generically for any editable field (z, t, T etc).
    
    % If val_flag is true
    if val_flag
        % Make a temp array of all values
        tick_temp = 1:upper_lim;

        % Find the discrete value closest to the entered number
        [~, minIdx] = min(abs(event.Value - tick_temp(:)));

        % Override the selected value
        x = tick_temp(minIdx);
    else
        % determine which discrete option the current value is closest
        % to when a value is entered, using the major tick
        [~, minIdx] = min(abs(event.Value - sld_x.MajorTicks(:)));

        % Override the selected value
        x = sld_x.MajorTicks(minIdx);
    end
    
    % Check the current position against the limits to determine if plus
    % and minus buttons are clickable
    check_limits(x, upper_lim, plus_x, minus_x);

    % Set the value of the slider to the current value
    sld_x.Value = x;

    % If an editable field exists
    if ~isempty(edt_x)
        % Set the value in the editable field to the currently selected
        % value
        edt_x.Value = x;
    end
end

function x = plusminusXButtonPushed(x, upper_lim, plus_x, minus_x,...
    sld_x, edt_x, plus_or_minus)
%PLUSMINUSXBUTTONPUSHED Updates x if a < or > button is clicked
%   
%   Input
%       x: value being changed
%       upper_lim: maximum or upper limit for x value
%       plus_x: uibutton handle for > button for x
%       minus_x: uibutton handle for < button for x
%       sld_x: slider handle for x
%       edt_x: editable field handle for x
%       plus_or_minus: 'plus' or 'minus' for increasing or decreasing by 1
% 
%   Output
%       x: updated value for x
%
%   Overview
%       Updates the value of x when a plus or minus button is clicked. In
%       addition, the function makes sure the  plus and minus buttons are
%       available depending on if the max or min is reached. The slider is
%       also updated. This function works generically for any value (z, t,
%       T etc). Disabling the button is included twice in case the user
%       clicks too fast.

    % Check the current position against the limits to determine if plus
    % and minus buttons are clickable
    check_limits(x, upper_lim, plus_x, minus_x);
    
    % If > was clicked
    if strcmp(plus_or_minus, 'plus')
        % Increase x by one
        x = x + 1;
    % Else if < was clicked
    elseif strcmp(plus_or_minus, 'minus')
        % Decrease x by 1
        x = x - 1;
    end
    
    % Check the current position against the limits to determine if plus
    % and minus buttons are clickable
    check_limits(x, upper_lim, plus_x, minus_x);
    
    % Set the slider value to the new value of x
    sld_x.Value = x;

    % If an editable field exists
    if ~isempty(edt_x)
        % Set the value in the editable field to the new value of x
        edt_x.Value = x;
    end
end

function check_limits(x, upper_lim, plus_x, minus_x)
%CHECK_LIMITS Checks to see if limits of x are reached and disables buttons
%   
%   Input
%       x: value being changed
%       upper_lim: maximum or upper limit for x value
%       plus_x: uibutton handle for > button for x
%       minus_x: uibutton handle for < button for x
% 
%   Output
%       None
%
%   Overview
%       Checks if the new value of x reaches the limits of x, and disables
%       the corresponding < or > button if the maximum or minimum is
%       reached.

    % If x is less than the max and greater than or equal to 1 (the
    % minimum)
    if (x < upper_lim) && (x >= 1)
        % Set the > button enable to true so it can be clicked
        plus_x.Enable = true;
    else
        % Set > button enable to false so it cannot be clicked
        plus_x.Enable = false;
    end

    % If x is greater than 1 (the minimum) and less than or equal to the
    % max
    if (x > 1) && (x <= upper_lim)
        % Set the < button enable to true so it can be clicked
        minus_x.Enable = true;
    else
        % Set the < button enbable to false so it cannot be clicked
        minus_x.Enable = false;
    end
end

function x = change_exp(sld_x, upper_lim, plus_x, minus_x,...
    edt_x, tick_flag)
%CHANGE_EXP Resets the sliders and buttons when the experiment is changed
%   
%   Input
%       sld_x: the slider for x
%       upper_lim: maximum or upper limit for x value
%       plus_x: uibutton handle for > button for x
%       minus_x: uibutton handle for < button for x
%       edt_x: entry field for editing the slider x
%       tick_flag: true false for setting the tick marks on the slider
% 
%   Output
%       x: The value being reset, will always be 1
%
%   Overview
%       This function will reset the slider and associated buttons/fields
%       when the experiment is changed. Each slider must be reset, except
%       brightness.
    
    % Set the limits for the slider to those of the new experiment
    sld_x.Limits = [1, upper_lim];
    edt_x.Limits = [1, upper_lim];
    
    % If true to set tick marks
    if tick_flag
        % Major and minor tick marks are set so that the labels are spread
        % out and rounded to the nearest 25
        tMinTick = ceil(upper_lim/350) * 5;
        tMajTick = tMinTick * 5;
        sld_x.MinorTicks = [1, tMinTick:tMinTick:upper_lim];
        sld_x.MajorTicks = [1, tMajTick:tMajTick:upper_lim];
    % Else, just set the tick mark step size to be 1
    else
        sld_x.MajorTicks = 1:upper_lim;
    end

    % Initilize x to one, set the slider to this value and set the plus
    % button so it is clickable and the minus button so it isn't
    % clickable
    x = 1;
    sld_x.Value = x;
    plus_x.Enable = true;
    minus_x.Enable = false;

    % If an editable field exists
    if ~isempty(edt_x)
        % Set the value in the editable field to the new value of x
        edt_x.Value = x;
    end
end