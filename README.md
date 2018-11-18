# imvol
Matlab-based Image Stack Viewer. Automatic ROI segmentation using local thresholding algorithm with a flexible update. Developed for analyzing time-varying fluorescence imaging data of a neural activity.


    Display using imshow() with interactive keyboard navigation for volume images (or stack)
    Select ROIs using imbinarize() with adjustable parameters and keyboard interactions. 
    (New figure will be created unless fig or axes handles are given.)

    Input: 
           vol - 2-D or 3-D matrix (stack)

    Varargin inputs: (optional)
          'title' - name string for experiment
          'roi'    - predefined cc structure for ROI mode. Final ROI will
          be differnet depending on your Sensitivity and Connectivity
          values. If you want to just repeat the input cc, add 'edit'
          option with 'false'. 

              Example:
                      imvol(vol, 'roi', cc, 'edit', false)
              
          'edit'   - When it is flase, the BW will not be modified by the mask or added ROIs manually drawn.
          'z_step_um' - z-stack spatial spacing between frames. Default is 1 um.
          'FOV'    - size of the image in um. Scale bar can be displayed.

    Output:
          hfig - ROI mask data will be saved in UserData field ((hfig.UserData.cc) 
          as well as in WorkSpace
                  

    Press 'spacebar' to switch between Modes.

       Mode1 - Imaging mode (default).
                  L /R    arrow keys - previous/next frames in vol images
                  Up/Down arrow keys - adjust contrast levels
                  'l' key - Draw line and get Z (or vertical) profile.
                            (default z_step = 1 um)
                  'b' key - Display scale bar. (Currently for 25x Leica obj).
                  's' key - Save current image as PNG
                  'v' key - Display/Hide verbose text notes in imshow
                  'g' key - Save image stack as GIF, MP4 and Tif formats.
                  'p' key - Create projected image. Max and Mean.

       Mode2 - ROI select mode. 'cc' (ROI) strucrue will be saved in Workspace.
                  L /R    arrow keys - adjust 'sensitivity' in imbinarize(J, 'adaptive')
                  Up/Down arrow keys - adjust 'connectivity threshold' in bwareaopen(bw, P_connected)

                  [keys for modifying ROIs]

                  'a' key - add whtie pixels (ellipse) in order to be detected as ROI.
                  'd' key - delete ROIs in drawed squre.
                  'l' key - add black pixels along the line to split ROIs.

                  [keys for visualizing ROIs]

                  'c' key - false color visualization for ROIs

    (c) 2018 Juyoung Kim 
