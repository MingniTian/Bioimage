% POIRegistration_103.m

% This script applies a geometric transformation on a selected image to
% register it to a selected reference image. The type of transformation is
% selected by the user as one of the following: 'Non-reflective
% similarity', 'Similarity', 'Affine' or 'Projective'. The user needs to
% manually pick several pairs of points of interest on the image to be
% transformed and on the reference image and then the program will do the
% transformation based on the spatial co-ordinates of the point pairs. The
% transformed image is saved as 'Transformed.tif', the merged image is
% saved as 'Merged.tif' and the fitting results are save in 'Fitting.txt'.

% Written in 2019 by Buyun Tian @ Tao Xu's lab, Institute of Biophysics,
% Beijing, China.

% v1.01 (2019) updated:
% Displayed image size and line colour optimised when selecting points
% (requires ginput_orange.m, a modified version of ginput.m)
% Allows user to adjust parameters for image merging until satisfied

% v1.02 (2020) updated:
% Outputs the transformation parameters (Transformation.mat) for later use
% Allows loading pre-defined transformation parameters when POI selection
% is not desired

% v1.03 (2020) updated:
% Allows zooming before point selection
% Bug fixed: if the image to be transformed is greyscale and the
% reference image is RGB or the other way round, there is no longer an
% error when the merged image is generated

    clear;
    close all;
    clc;
    
    % Select transformation type
    list = {'Non-reflective similarity', 'Similarity', 'Affine', 'Projective'};
    [indx,tf] = listdlg('ListString', list, 'PromptString', 'Select a transformation type', 'SelectionMode', 'single', 'ListSize', [260,120], 'InitialValue', 4, 'Name', 'Transformation types');
    if tf == 0
        uiwait(errordlg('No transformation type selected.', 'Error'));
        return;
    end
    
    % Read image files
    uiwait(helpdlg('Select image files.','Instruction'));
    [filename1, filepath1] = uigetfile('*.*','Select image to be transformed');
    OrigImg = im2double(imread([filepath1 filename1])); % Original image to be transformed
    [filename2, filepath2] = uigetfile('*.*','Select reference image');
    RefImg = im2double(imread([filepath2 filename2]));  % Reference image
    
    % Display images
    ScreenSize = get(0,'ScreenSize');
    ScreenWidth = ScreenSize(3);
    ScreenHeight = ScreenSize(4);
    OrigImgSize = size(OrigImg);
    RefImgSize = size(RefImg);  
    OrigImgMag = 80*min(ScreenWidth/OrigImgSize(2), ScreenHeight/OrigImgSize(1));  % Magnification of original image when displayed
    RefImgMag = 80*min(ScreenWidth/RefImgSize(2), ScreenHeight/RefImgSize(1));  % Magnification of reference image when displayed
    figure, imshow(OrigImg, [], 'InitialMagnification', OrigImgMag);
    title('Image to be transformed');
    figure, imshow(RefImg, [], 'InitialMagnification', RefImgMag);
    title('Reference image');
    
    % Choose whether to load a pre-defined transformation
    answer = questdlg('Load a pre-defined transformation?', 'Processing mode', 'Yes', 'New transformation', 'New transformation');
    switch answer
        
    % For pre-defined transformation
        case 'Yes'
            [TrfFileName, TrfFilePath] = uigetfile('*.mat','Select transformation file');
            load([TrfFilePath TrfFileName]);
            
    % For POI transformation
        case 'New transformation'
            
            % Manually select matching points
            uiwait(helpdlg('Press any key to start selection once zooming is OK. Click to select points of interest. Press Enter when you finish.','Instruction'));
            POI_O = [];	% Points of interest on original image
            POI_R = [];	% Points of interest on reference image
            for i = 1:256
                figure(1);
                hold on;
                zoom on;
                pause;
                zoom off;
                [x1, y1] = ginput_orange(1);	% Select a point on original image
                zoom out;
                if (isempty(x1) || isempty(y1))
                   break;
                end
                scatter(x1, y1, 'r', '.');  % Mark selected point
                text(x1, y1, [' ', num2str(i)], 'Color', 'red');
                hold off;
                figure(2);
                hold on;
                zoom on;
                pause;
                zoom off;
                [x2, y2] = ginput_orange(1);   % Select a corresponding point on reference image
                zoom out;
                if (isempty(x2) || isempty(y2))
                    break;
                end
                scatter(x2, y2, 'r', '.');  % Mark selected point
                text(x2, y2, [' ', num2str(i)], 'Color', 'red');    
                hold off;
                POI_O = [POI_O; x1 y1];
                POI_R = [POI_R; x2 y2];
            end
    
            % Error checkpoint 
            N = size(POI_O,1);  % Number of points of interest on original image
            if N ~= size(POI_R,1)
                uiwait(errordlg('Points do not match.', 'Error'));
                return;
            end
            N_min = 0;  % Minimum number of pairs of points of interest
            TrfType = 'None';   % Transformation type
            switch indx
                case 1
                    TrfType = 'nonreflectivesimilarity';
                    N_min = 2;
                case 2
                    TrfType = 'affine';
                    N_min = 3;
                case 3
                    TrfType = 'similarity';
                    N_min = 3;
                case 4
                    TrfType = 'projective';
                    N_min = 4;
                otherwise
                    uiwait(errordlg('Invalid transformation type.', 'Error'));
                    return;
            end 
            if N < N_min
                uiwait(errordlg('Not enough points.', 'Error'));
                return;
            end
    
            % Compute & save transformation
            tform = fitgeotrans(POI_O, POI_R, TrfType);
            save('Transformation.mat', 'tform');
    
            % Check if points of interest are translocated properly
            [tf_pt_x, tf_pt_y] = transformPointsForward(tform, POI_O(:, 1),POI_O(:, 2));    % Transformed x, y co-ordinates of points of interest
            figure(2);
            hold on;
            scatter(tf_pt_x, tf_pt_y, 'g', '.');    % Mark transformed positions of points of interest
            hold off;
            helpdlg('Transformed positions of points of interest are marked as green dots on the reference image.','Message')
    end
        
    % Transform original image and show results
    TrfImg = imwarp(OrigImg, tform, 'cubic', 'OutputView', imref2d(size(RefImg)));  % Transformed image
    figure, imshow(TrfImg, []);
    title('Transformed image');
    figure, imshowpair(RefImg, TrfImg);
    title('Comparison');
    figure;
    title('Merged image');
    MergImg = [];  % Merged image
    while 1  % Allows user to adjust merging parameters until the process is manually cancelled
        DispPara = inputdlg({'Weight of transformed image:', 'Weight of reference image:'}, 'Adjust parameters for merged image until cancelled', 1, {'1', '1'}, 'on');
        try
            weight1 = str2double(DispPara{1});
            weight2 = str2double(DispPara{2});
            if (numel(size(TrfImg)) == 3) && (numel(size(RefImg)) == 2)
                MergImg = imlincomb(weight1, TrfImg, weight2, cat(3, RefImg, RefImg, RefImg));
            end
            if (numel(size(TrfImg)) == 2) && (numel(size(RefImg)) == 3)
                MergImg = imlincomb(weight1, cat(3, TrfImg, TrfImg, TrfImg), weight2, RefImg);
            end
            if numel(size(TrfImg)) == numel(size(RefImg))
                MergImg = imlincomb(weight1, TrfImg, weight2, RefImg);
            end
            imshow(MergImg, []);
        catch
            break;
        end
    end   
    
    % Save images
    imwrite(im2uint16(TrfImg), 'Transformed.tif');
    imwrite(im2uint16(MergImg), 'Merged.tif');

    % Evaluate fitting
    switch answer
        case 'Yes'
            % No fitting
        case 'New transformation'
            % Compute parameters for evaluation
            POI_T = [tf_pt_x tf_pt_y];  % Transformed positions of points of interest
            Centre = mean(POI_R);   % Variable used to compute R-squared
            Residuals = POI_R - POI_T;
            ResidualMag = zeros(N,1);   % Magnitudes of residuals
            RelativeLoc = POI_R;    % Variable used to compute R-squared
            for i = 1:N
                ResidualMag(i) = norm(Residuals(i,:));
                RelativeLoc(i,:) = RelativeLoc(i,:) - Centre;
            end
            SSE = sum(ResidualMag.^2);  % Sum of squares due to error
            MSE = SSE/N;    % Mean squared error
            RMSE = MSE^0.5; % Root mean squared error
            SST = sum(sum(RelativeLoc.^2)); % Total sum of squares
            R_squared = 1 - SSE/SST;
    
            % Save evaluation of fitting
            fid = fopen('Fitting.txt', 'w');
            fprintf(fid, '%s\r\n', ['Transformation type: ' lower(list{indx})]);
            fprintf(fid, '%s\r\n', ['Number of pairs of points of interest: ' num2str(N)]);
            fprintf(fid, '%s\r\n', ['SSE, sum of squares due to error: ' num2str(SSE)]);
            fprintf(fid, '%s\r\n', ['MSE, mean squared error: ' num2str(MSE)]);
            fprintf(fid, '%s\r\n', ['RMSE, root mean squared error: ' num2str(RMSE)]);
            fprintf(fid, '%s\r\n', ['R_squared: ' num2str(R_squared)]);
            fprintf(fid, ['Magnitudes of residuals: ' num2str(ResidualMag')]);
            fclose(fid);
    end
