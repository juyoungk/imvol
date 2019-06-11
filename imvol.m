function [hfig] = imvol(vol, varargin)
%IMVOL Display volume images and select ROIs
%     Display using imshow() with interactive keyboard navigation for volume images (or stack)
%     Select ROIs using imbinarize() with adjustable parameters and keyboard interactions. 
%     (New figure will be created unless fig or axes handles are given.)
% 
%     Input: 
%            vol - 2-D or 3-D matrix (stack)
% 
%     Varargin inputs: (optional)
%           'title' - name string for experiment
%           'roi'    - predefined cc structure for ROI mode. Final ROI will
%           be differnet depending on your Sensitivity and Connectivity
%           values. If you want to just repeat the input cc, add 'edit'
%           option with 'false'. 
% 
%               Example:
%                       imvol(vol, 'roi', cc, 'edit', false)
%               
%           'edit'   - When it is flase, the BW will not be modified by the mask or added ROIs manually drawn.
%           'z_step_um' - z-stack spatial spacing between frames. Default is 1 um.
%           'FOV'    - size of the image in um. Scale bar can be displayed.
%
% 
%     Output:
%           hfig - ROI mask data will be saved in UserData field ((hfig.UserData.cc) 
%           as well as in WorkSpace.
%                   
% 
%     Press 'spacebar' to switch between Modes.
% 
%        Mode1 - Imaging mode (default).
%                   L /R    arrow keys - previous/next frames in vol images
%                   Up/Down arrow keys - adjust contrast levels
%                   'l' key - Draw line and get Z (or vertical) profile.
%                             (default z_step = 1 um)
%                   'r' key - Draw rectangular and get Z (or vertical) profile of Max projection line.
%                             (default z_step = 1 um)
%                   'b' key - Display scale bar. (Currently for 25x Leica obj).
%                   's' key - Save current image as PNG
%                   'v' key - Display/Hide verbose text notes in imshow
%                   'g' key - Create animated movies as GIF, MP4 and Tif
%                   formats.
%                   'p' key - Create projected image. Max and Mean.
%                   'z' key - display Z relative position of current frame (Default z-step is 1 um).                  
% 
%        Mode2 - ROI select mode. 'cc' (ROI) strucrue will be saved in Workspace.
%                   L /R    arrow keys - adjust 'sensitivity' in imbinarize(J, 'adaptive')
%                   Up/Down arrow keys - adjust 'connectivity threshold' in bwareaopen(bw, P_connected)
% 
%                   [keys for modifying ROIs]
% 
%                   'a' key - add whtie pixels (ellipse) in order to be detected as ROI.
%                   'd' key - delete ROIs in drawed squre.
%                   'l' key - add black pixels along the line to split ROIs.
% 
%                   [keys for visualizing ROIs]
% 
%                   'c' key - false color visualization for ROIs
% 
%     (c) 2018 Juyoung Kim 

    p = ParseInput(varargin{:});
    %
    s_title = p.Results.title;
    
    hfig = p.Results.hfig;
    ax = p.Results.axes;
    cc = p.Results.roi;
    zoom = p.Results.scanZoom;
    SAVE_png = p.Results.png;
    FLAG_txt = p.Results.verbose;
    FLAG_z = p.Results.z_display;
    FLAG_scale_bar = true;
    FLAG_roi = false;
    FLAG_color_segmentation = false;
    FLAG_hole_fill = true;
    FLAG_edit = p.Results.edit;
    globalContrast = p.Results.globalContrast;
    z_step_um = p.Results.z_step_um;
    timestamp = p.Results.timestamp;
    %
%     import java.awt.Robot
%     import java.awt.event.*
%     keys = Robot;
%     keys.setAutoDelay(0);
    %
    if nargin < 1
        error('No image (stack) was not provided.');
    else
        vol_inputname = inputname(1);
    end
    
    % str of input variable
    if isempty(s_title)
        s_title = vol_inputname;
    end
    % replace '_' with ' '
    s_title = strrep(s_title, '_', ' ');
    
    N = ndims(vol);
    if N > 3
        error('The input image stack (vol) has too high dims >3.');
    elseif N < 2
        error('Not image (ndims <2)');
    end
    
    if ishandle(hfig)
        % if fig handle is given, give focus to the figure.
        figure(hfig);
    elseif ishandle(ax)
        % if axes handle is given, figure out the container handle.
        axes(ax);
        hfig = ax.Parent;
    else
        if ~isempty(hfig)
            disp('(imvol) The input fig handle was not appropriate. New figure was created.');
        end
        hfig = figure();
    end
    hfig.Color = 'none';
    hfig.PaperPositionMode = 'auto';
    hfig.InvertHardcopy = 'off';   
    pos = hfig.Position;
    axes('Position', [0  0  1  0.9524], 'Visible', 'off');
    
    % ex str in UserData
    %hfig.UserData.ex_str = ex_str;
    
    % Set the callback on figure
    set(hfig, 'KeyPressFcn', @keypress)
    
    % Get frame numbers
    [rows, cols, n_frames] = size(vol);
    
    % Aux info: timestamp
    if isempty(timestamp)
        timestamp = NaN(1, n_frames);
    elseif n_frames ~= length(timestamp)
        disp('num of timestamp ~= num of frames.');
        timestamp = NaN(1, n_frames);
    end
    
    % mask variables for ROI removal by user
    mask = false(rows, cols);
    white = false(rows, cols);
    r = []; c = []; % row, col for selected points by clicks
    
    % Default parameters
    data.i = 1; % index for stack
    data.imax = n_frames;
    tols = [0, 0.01, 0.02, 0.05, 0.1:0.1:0.9, 1:0.2:2, 2.5:0.5:5, 6:1:11, 12:2:20, 25:5:95]; % percentage; tolerance for saturation
    n_tols = length(tols);
    id_tol = 4; % initial tol = 0.05;
    id_add_lower = 1;
    id_add_upper = 3;
    
    % ROI mode parameters
    sensitivity_0 = 0.04; % sensitivity for adaptive binarization
    P_connected_0 = 55; % depending on magnification (zoom) factor
    sensitivity = sensitivity_0; 
    P_connected = P_connected_0; 
    
    % Global MinMax
    if globalContrast
        vv = vol(:);
        minI = min(vv);
        maxI = max(vv);
        vv = mat2gray(vv);
    end
        
    % Nested function definition for easy access to stack 'vol'
    function redraw()
        % get focus
        figure(hfig);
        if N == 2
            I = vol;
        else
            I = comp(vol, data.i);
        end
        
        % Contrast range
        upper = max(1 - (tols(id_tol) + tols(id_add_upper))*0.01, 0);
        lower = min((tols(id_tol) + tols(id_add_lower))*0.01, upper);
        Tol = [lower upper];
        
        % MinMax
        if globalContrast
            I = mat2gray(I);
            MinMax = stretchlim(vv, Tol); % Can vary depending on Tol values. 
        else
            minI = min(I(:));
            maxI = max(I(:));
            % Rescaling to [0 1] for stretchlim
            I = mat2gray(I);
            MinMax = stretchlim(I,Tol); % output is always [0 1].
        end
        % Image rescaling
        J = imadjust(I, MinMax);
        
        % Reconversion to original scale (for display)
        MinMax = MinMax * double(maxI - minI);
        MinMax = MinMax + double(minI);
        
        % new way of adjusting?
        %J = imadjust(I, MinMax, [0 1]); % MinMax can be real value?

        % draw image or visualization
        if ~FLAG_roi 
            imshow(J);
            
            % turn on this option for 2D presentation of 1D bar patterns
            % daspect auto;
            
            % txt str
            %str1 = sprintf('low=%.3f upp=%.3f', lower, upper);
            str1 = sprintf('Min =%6.0f (%4.1f%%)\nMax =%5.0f (%4.1f%%)',MinMax(1),lower*100,MinMax(2),upper*100);
            %str2 = '''q'' for default contrast. ''SPACE'' for ROI mode. ''b'' scale bar';
            str2 = '''SPACE'' for ROI mode.';
            str3 = sprintf('%d/%d', data.i, data.imax);
        else 
            % ROI mode
            if ~isempty(p.Results.roi)
                % if cc is given
                cc = p.Results.roi;
                bw = conn_to_bwmask(cc);
                bw = max(bw, [], 3);
            else
                bw = imbinarize(J, 'adaptive', 'Sensitivity', sensitivity);
            end
            
            if FLAG_edit
            % edit ROI. Update bw
                bw = bw & (~mask);    % get ROI mask and then subtract it from the image
                bw = bw | (white);
                bw = bwareaopen(bw, P_connected); % remove small area
                if FLAG_hole_fill
                    bw = imfill(bw, 'hole');
                end
                bw = bw - bwselect(bw, c, r, 8);  % remove mouse-clicked components
                %s = regionprops(bw, 'Centroid');
            end
            cc = bwconncomp(bw, 8); % 'cc' is updated inside a local function.
            
            % save for regenrating same pattern
                cc.mask  = mask;
                cc.white = white;
                cc.P_connected = P_connected;
                cc.sensitivity = sensitivity;
            % save frame id & aux info
                cc.i_image = data.i;
                cc.timestamp = timestamp(data.i);
                cc.numImages = data.imax;

            % visualization of computed ROI
            if ~FLAG_color_segmentation
                imshow(J); 
                
                % n_frames =2 case: 
                % imshowpair?
                
                hold on
                    % Contour 
                    visboundaries(bw,'Color','r','LineWidth', 0.7); 

                    % ROI number display
                    s = regionprops(cc, 'extrema');
                    for k = 1:numel(s)
                       e = s(k).Extrema;
                       text(e(4,1), e(4,2), sprintf('%d', k), 'Color', 'r', ... %5th comp: 'bottom-right'
                           'VerticalAlignment', 'bottom', 'HorizontalAlignment','left'); 
                    end
                hold off
            else
                labeled = labelmatrix(cc);
                RGB_label = label2rgb(labeled, @parula, 'k', 'shuffle');
                imshow(RGB_label);
            end
            % update 'cc' whenever ROI mode is activated.
                % figure handle
                hfig.UserData.cc = cc;
                % in workspace
                v_name = 'cc';
                assignin('base', v_name, cc);
                assignin('base', 'mask', mask);
                assignin('base', 'white', white);
                %
                disp([num2str(cc.NumObjects), ' Objects (ROIs) are selected.']);    
        end
        
        % Shared components
        ax = gca;
        % Title (on axes, not image)
        title(s_title, 'FontSize', 15, 'Color', 'w');
        
        % Scale bar
        if FLAG_scale_bar
            fov = p.Results.FOV;
            if fov > 0 
                % Positive if FOV is given directly as varargin by user.
                % do nothing
            else
                % estimation by zoom factor
                fov = get_FOV_size_x25_Leica(zoom);
            end
            
            if fov > 0 % Draws scale bar only when fov > 0
                px_per_um = rows/fov;
                hold on;
                    if zoom < 2
                        l_scalebar = 100; % um
                        x0 = ax.XLim(end) * 0.75;
                        y0 = ax.YLim(end) * 0.94;
                    else
                        l_scalebar = 30; % um
                        x0 = ax.XLim(end) * 0.80;
                        y0 = ax.YLim(end) * 0.94;
                    end
                    quiver(x0, y0, l_scalebar*px_per_um, 0, 'ShowArrowHead', 'off', 'Color', 'w', 'LineWidth', 2);
                    text(x0+l_scalebar*px_per_um/2, y0, [num2str(l_scalebar),' um'], 'FontSize', 15, 'Color', 'w', ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment','center');
                hold off;
                %line( [100 200], round(0.85*ax.YLim(end)), 'Color', 'w', 'LineWidth', 4);
            end
        end
        
        if SAVE_png
            filename = strrep(s_title, ' ', '_');
            filename = strrep(filename, '(', '_');
            filename = strrep(filename, ':', '');
            if ~FLAG_roi
                filename = sprintf('%s_%dof%d_lower%.3f_upper%.3f.png', filename, data.i, n_frames, lower, upper); 
            else
                filename = sprintf('%s_%dof%d_ROI_conn%.3f_sens%.3f.png', filename, data.i, n_frames, P_connected, sensitivity);
            end
            %saveas(hfig, filename);
            print(filename, '-dpng', '-r300'); %high res
            SAVE_png = false; % save only one time
        end   
        
        if FLAG_txt
            % Text on the image. (will be cleared if redrawed.)
            if FLAG_roi
                str1 = sprintf('Sens.=%.2f Pconn.=%.2f.', sensitivity, P_connected);
                str2 = sprintf('Press ''d'' or ''r'' to remove ROIs. ''q'' for default settings');
                str3 = sprintf('%d/%d', data.i, data.imax);
%             else
%                 str1 = sprintf('low=%.3f upp=%.3f', lower, upper);
%                 str2 = '''q'' for default contrast. ''SPACE'' for ROI mode. ''b'' scale bar';
%                 str3 = sprintf('%d/%d', data.i, data.imax);
            end
            
            % x,y for text. Coordinate for imshow is different from plot
            text(ax.XLim(1), ax.YLim(end), str1, 'FontSize', 14, 'Color', 'w', ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment','left');
            text(ax.XLim(2), ax.YLim(end), str3, 'FontSize', 14, 'Color', 'w', ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment','right');
            text((ax.XLim(1)+ax.XLim(end))/2, ax.YLim(end), str2, 'FontSize', 14, 'Color', 'w', ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment','center');
        end
        
        if FLAG_z
            %str3 = sprintf(' (Frame %d/%d)\n  Z  = %3d um ', data.i, data.imax, data.i*z_step_um);
            str3 = sprintf('   Z =%3d um ', data.i*z_step_um);
            text(ax.XLim(1), ax.YLim(1), str3, 'FontSize', 15, 'Color', 'w', ...
                'VerticalAlignment', 'top', 'HorizontalAlignment','left');
        end
        
        uiresume(hfig);
    end

    % Update the display of the surface
    redraw();
    
    % Normal mode: viewer
    function keypress(~, evnt)
        
        switch lower(evnt.Key)
            case 'rightarrow'
                data.i = min(data.i + 1, data.imax); 
            case 'leftarrow'
                data.i = max(1, data.i - 1);
            case 'uparrow'
                id_tol = min(id_tol + 1, n_tols);
            case 'downarrow'
                id_tol = max(1, id_tol - 1); 
            case '1'
                id_add_lower = max(1, id_add_lower - 1); 
            case '2'
                id_add_lower = min(id_add_lower + 1, n_tols);
            case '9'
                id_add_upper = max(1, id_add_upper - 1); 
            case '0'
                id_add_upper = min(id_add_upper + 1, n_tols);
            case 'l' % line profile
                [~,~,c,xi,yi] = improfile;
                c_section = zeros(length(c), n_frames);
                for k = 1:n_frames
                    c_section(:,k) = improfile(vol(:,:,k),xi,yi);
                end
                if zoom == 0
                    zoom = input('Please enter the imaging zoom factor: ');
                end
                fov = get_FOV_size_x25_Leica(zoom);
                px_per_um = rows/fov;
                a_ratio = px_per_um/z_step_um;
                img = c_section.';
                [numrows, numcols] = size(img);
                C = imresize(img, [round(a_ratio*numrows), numcols]);
                make_im_figure(500, 0); 
                myshow(C, 0.1);                
                ax = gca; 
                %scale bar?l
                hold on;
                l_scalebar = 50; % um
                x0 = ax.XLim(end) * 0.90;
                y0 = ax.YLim(end) * 0.90;
                quiver(x0, y0, l_scalebar*px_per_um, 0, 'ShowArrowHead', 'off', 'Color', 'w', 'LineWidth', 2);
                text(x0+l_scalebar*px_per_um/2, y0, [num2str(l_scalebar),' um'], 'FontSize', 15, 'Color', 'w', ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment','center');
                hold off;
            case 'r' % rectangular profile
                h = imrect;
                m = h.createMask;
                [row, col] = size(m);
                c_section = zeros(col, n_frames);
                for k = 1:n_frames
                    im = vol(:,:,k);
                    im = im .* int16(m);r
                    %im = myadjust(im, 0.05);
                    %im = im .* double(m);
                    %c_section(:,k) = mean(im, 1);
                    c_section(:,k) = max(im, [], 1);
                end
                if zoom == 0
                    zoom = input('Please enter the imaging zoom factor: ');
                end
                fov = get_FOV_size_x25_Leica(zoom);
                px_per_um = rows/fov;
                a_ratio = px_per_um/z_step_um;
                img = c_section.';
                [numrows, numcols] = size(img);
                C = imresize(img, [round(a_ratio*numrows), numcols]);
                make_im_figure(500, 0); 
                myshow(C, 0.3);                
                ax = gca; 
                %scale bar?
                hold on;
                l_scalebar = 50; % um
                x0 = ax.XLim(end) * 0.80;
                y0 = ax.YLim(end) * 0.90;
                quiver(x0, y0, l_scalebar*px_per_um, 0, 'ShowArrowHead', 'off', 'Color', 'w', 'LineWidth', 2);
                text(x0+l_scalebar*px_per_um/2, y0, [num2str(l_scalebar),' um'], 'FontSize', 15, 'Color', 'w', ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment','center');
                hold off;
            case 's' % save current snapshot
                SAVE_png = true;
            case 'g' %save as GIF and MP4
                % conditions
                FLAG_txt = false;
                if contains(s_title, 'stack')
                    FLAG_z = true;
                else
                    FLAG_z = false;
                end
                delaytime = 0.1;
                %
                v = VideoWriter([s_title, '.mp4'], 'MPEG-4'); 
                v.FrameRate = 4;
                open(v);
                % 
                for k = 1:data.imax
                    data.i = k;
                    redraw(); title('');                  
                    frame = getframe(hfig); % frame from Handle hfig.
                    
                    % animated GIF
                    im = frame2im(frame);
                    [A, map] = rgb2ind(im, 256);
                    if k == 1 
                        imwrite(A, map, [s_title, '.gif'], 'gif', 'LoopCount', Inf, 'DelayTime', delaytime);
                        %imwrite(im, [s_title, '.tif']);
                    else
                        imwrite(A, map, [s_title, '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', delaytime);
                        %imwrite(im, [s_title, '.tif'], 'WriteMode', 'append');
                    end
                    % 
                    writeVideo(v, frame);
                end
                close(v);
                FLAG_txt = true;
                
            case 'p' % projection
                prompt = {'Start layer: ','Last layer: ','Projection type: '};
                d_title = 'Projection image';
                dims = [1 15];
                definput = { num2str(data.i), num2str(data.i), 'Max'};
                answer = inputdlg(prompt,d_title,dims,definput);
                % check if out-of-range?
                img = vol(:,:,str2double(answer{1}):str2double(answer{2}));
                hnew = make_im_figure(pos(3), 0);
                imvol( max(img, [], 3), 'title', [answer{3}, ' projection (', answer{1},'-',answer{2},')'], 'hfig', hnew, 'scanZoom', zoom);
                hnew = make_im_figure(2*pos(3)-pos(1), 0);
                imvol( mean(img, 3), 'title', ['mean', ' projection (', answer{1},'-',answer{2},')'], 'hfig', hnew, 'scanZoom', zoom);
                
            case 'b'
                FLAG_scale_bar = ~FLAG_scale_bar;
                % ask FOV size if there no function for get_FOV_size
                %
            case 'space' % ROI mode switch
                FLAG_roi = ~FLAG_roi;
                set(hfig, 'KeyPressFcn', @keypress_roi)
            case 'v' % verbose output
                FLAG_txt = ~FLAG_txt;
            case 'q' % default contrast
                id_tol = 5;
                id_add_lower = 1;
            case 'z' % z depth display mode
                FLAG_z = ~FLAG_z;
                if FLAG_z
                    % display conversion unit
                    fprintf('Z step between frames is %.2f um.\n', z_step_um);
                end
            otherwise
                uiresume(hfig)
                return;
        end
        uiresume(hfig)
        redraw();
    end
    
    % ROI mode
    function keypress_roi(~, evnt)
        % step for tolerance or contrast
        s = 0.02;
        s_pixel = 5;
        
        % ESC key press? does not work for imroi..
%         keys.keyPress(java.awt.event.KeyEvent.VK_ESCAPE) 
%         keys.keyRelease(java.awt.event.KeyEvent.VK_ESCAPE)
%         pause(0.05)
  
        switch lower(evnt.Key)
            case 'rightarrow'
                sensitivity = min(sensitivity + s, 1); 
            case 'leftarrow'
                sensitivity = max(0, sensitivity - s); 
            case 'uparrow' % ignore disconnected dots. remove noise.
                P_connected = P_connected + s_pixel;
            case 'downarrow'
                P_connected = max(0, P_connected - s_pixel);
            case '1'
                id_add_lower = max(1, id_add_lower - 1); 
            case '2'
                id_add_lower = min(id_add_lower + 1, n_tols);
            case '9'
                id_add_upper = max(1, id_add_upper - 1); 
            case '0'
                id_add_upper = min(id_add_upper + 1, n_tols);
            case 's'
                SAVE_png = true;
            case 'c' % Color map
                FLAG_color_segmentation = ~FLAG_color_segmentation;
            case 'space' 
                % Go back to imshow mode
                FLAG_roi = ~FLAG_roi;
                set(hfig, 'KeyPressFcn', @keypress)
            case 'r' % mask update. remove connected components by multiple mouse clicks
                uiresume(hfig)
                [col, row] = getpts;
                c = [c; col];
                r = [r; row];
                
            case 'd' % mask update. 'Drag': remove all components in specified rect ROI.
                %uiresume(hfig)
                hrect = imrect;
                while ~isempty(hrect)
                    m = createMask(hrect);
                    mask = mask | m;
                    redraw();
                    hrect = imrect;
                end
                
            case 'l' % line mask
                %uiresume(hfig)
                hline = imline;
                % additional callback to imline object? to prevent
                % crosstalk.. 
                m = createMask(hline);
                mask = mask | m;
                
%                 while ~isempty(hline)
%                     m = createMask(hline);
%                     mask = mask | m;
%                     redraw();
%                     hline = imline;
%                 end
  
            case 'a' % add patch
                uiresume(hfig)
                h_ellip = imellipse;
                m = createMask(h_ellip);
                white = white | m;         
            
            case 'f' % turn off automatic filling (imfill) inside the grain.
                FLAG_hole_fill = ~FLAG_hole_fill;
                
            case 'v' % verbose output
                FLAG_txt = ~FLAG_txt;
                
            case 'q' % default contrast
                sensitivity = sensitivity_0;
                P_connected = P_connected_0;
                mask = false(rows, cols);
                white = false(rows, cols);
                r = []; c = [];
           
            otherwise
                return; % must exit the function for other key presses.
        end
        
        uiresume(hfig)
        redraw();
    end

    
    
end

function q = is_valid_zoom(zoom)
    if nargin < 1
        q = false;
    end

    if (zoom >= 1) & (zoom <=7)
        q = true;
    else
        q = false;
    end
end

function hfig = make_im_figure(x_shift, y_shift)
    
    if nargin < 2
        y_shift = 0;
    end
    
    if nargin < 1
        x_shift = 0;
    end
    
    pos = get(0, 'DefaultFigurePosition');
    hfig = figure('Position', [pos(1) + x_shift, pos(2) + y_shift, pos(3), pos(4)]);
    
    hfig.Color = 'none';
    hfig.PaperPositionMode = 'auto';
    hfig.InvertHardcopy = 'off';   
    
    axes('Position', [0  0  1  0.9524], 'Visible', 'off'); % space for title
end

function [roi_array, bw_selected] = conn_to_bwmask(cc, id_selected)
% convert cc to bwmask array and totla bw for selected IDs.

    if nargin < 2
        id_selected = 1:cc.NumObjects;
    end

    roi_array = false([cc.ImageSize, cc.NumObjects]);

    for i = 1:cc.NumObjects
        grain = false(cc.ImageSize);
        grain(cc.PixelIdxList{i}) = true;
        roi_array(:,:,i) = grain;
    end
    % Total bw for selected IDs
    bw_selected = max( roi_array(:,:,id_selected), [], 3);
    
end

function p =  ParseInput(varargin)
    
    p  = inputParser;   % Create an instance of the inputParser class.
    
    addParamValue(p,'title', []);
    addParamValue(p,'hfig', []);
    addParamValue(p,'axes', []);
    addParamValue(p,'roi', []);
    %addParamValue(p,'sync', []);
    addParamValue(p,'verbose', true, @(x) islogical(x));
    addParamValue(p,'png', false, @(x) islogical(x));
    addParamValue(p,'scanZoom', 0, @(x) isnumeric(x));
    addParamValue(p,'edit', true, @(x) islogical(x));
    addParamValue(p,'z_step_um', 1, @(x) isnumeric(x));
    addParamValue(p,'z_display', false, @(x) islogical(x));
    addParamValue(p,'FOV', 0, @(x) isnumeric(x)); % um
    addParamValue(p,'timestamp', [], @(x) isnumeric(x)); % sec
    addParamValue(p,'globalContrast', false, @(x) islogical(x)); % um
    
    % Call the parse method of the object to read and validate each argument in the schema:
    p.parse(varargin{:});
    
end

function Xcomp = comp(X, range)
% Pick a specific range of frames or components of X
% , which is the last dim of matrix X (usually dim.3 for 2-D images)
% 
% if Dim ==2
%     Xcomp = X(:,range);
% elseif Dim ==3
%     Xcomp = X(:,:,range);
% elseif Dim == 4
%     Xcomp = X(:,:,:,range);
% elseif Dim == 5
%     Xcomp = X(:,:,:,:,range);
%

Dim = ndims(X);
d = size(X);

if max(range)>d(end) || min(range)<1
    disp('Function comp: components (frames) are out of range');
    Xcomp = [];
    return;
end

if Dim ==2 % 2-D matrix
    Xcomp = X(:,range);
elseif Dim ==3
    Xcomp = X(:,:,range);
elseif Dim == 4
    Xcomp = X(:,:,:,range);
elseif Dim == 5
    Xcomp = X(:,:,:,:,range);
else
    Xre = reshape(img,[],size(img,ndims(img)));
    Xcomp = Xre(:,range);
    disp('Dim of matrix X @ comp function is more than 5. X was reshaped into 1-D.');
end

end

function J = myshow(I, c)
%imshow with contrast value (%)
    J = myadjust(I, c);
    imshow(J);
end

function J = myadjust(I, c)
%imshow with contrast value (%)
    I = mat2gray(I);
    Tol = [c*0.01 1-c*0.01];
    MinMax = stretchlim(I,Tol);
    J = imadjust(I, MinMax);
end

% redefine imline? 
