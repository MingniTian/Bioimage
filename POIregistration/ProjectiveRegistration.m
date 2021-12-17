% ProjectiveRegistration.m

% This script applies a projective transformation on a selected image to
% register it to a selected reference image. Users need to manually pick
% several pairs of points of interest on the image to be transformed and on
% the reference image and then the program will do the transformation based
% on the spatial co-ordinates of the point pairs. The reprojected image is
% saved as 'Reprojected.tif', the merged image is saved as 'Merged.tif' and
% the fitting results are save in 'Fitting.txt'.

% Written in 2019 by Buyun Tian @ Tao Xu's lab, Institute of Biophysics,
% Beijing, China.

    clear;
    close all;
    clc;
    
    % Read image files
    uiwait(helpdlg('Select image files.','Instruction'));
    [filename1, filepath1]=uigetfile('*.*','Select image to be transformed');
    OrigImg = im2double(imread([filepath1 filename1])); % Original image to be transformed
    [filename2, filepath2]=uigetfile('*.*','Select reference image');
    RefImg = im2double(imread([filepath2 filename2]));  % Reference image
    
    % Manually select matching points
    figure, imshow(OrigImg, []);
    title('Image to be transformed');
    figure, imshow(RefImg, []);
    title('Reference image');
    uiwait(helpdlg('Click to select points of interest. Press Enter when you finish.','Instruction'));
    POI_O = [];	% Points of interest on original image
    POI_R = [];	% Points of interest on reference image
    for i = 1:256
        figure(1);
        hold on;
        [x1, y1] = ginput(1);	% Select a point on original image
        if (isempty(x1) || isempty(y1))
            break;
        end
        scatter(x1, y1, 'r', '.');  % Mark selected point
        text(x1, y1, [' ', num2str(i)], 'Color', 'red');
        hold off;
        figure(2);
        hold on;
        [x2, y2] = ginput(1);   % Select a corresponding point on reference image
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
    if N < 4
        uiwait(errordlg('Not enough points.', 'Error'));
        return;
    end
    
    % Compute transformation
    tform = fitgeotrans(POI_O, POI_R, 'projective');
    
    % Check if points of interest are translocated properly
    [tf_pt_x, tf_pt_y] = transformPointsForward(tform, POI_O(:, 1),POI_O(:, 2));    % Transformed x, y co-ordinates of points of interest
    figure(2);
    hold on;
    scatter(tf_pt_x, tf_pt_y, 'g', '.');    % Mark transformed positions of points of interest
    hold off;
    helpdlg('Transformed positions of points of interest are marked as green dots on the reference image.','Message')
    
    % Transform original image and show results
    TrfImg = imwarp(OrigImg, tform, 'cubic', 'OutputView', imref2d(size(RefImg)));
    figure, imshow(TrfImg, []);
    title('Transformed image');
    figure, imshowpair(RefImg, TrfImg);
    title('Comparison');
    DispPara = inputdlg({'Weight of reprojected image:', 'Weight of reference image:'}, 'Input parameters for merged image', 1, {'1', '1'}, 'on');
    weight1 = str2double(DispPara{1});
    weight2 = str2double(DispPara{2});
    MergImg = imlincomb(weight1, TrfImg, weight2, RefImg);
    figure, imshow(MergImg, []);
    title('Merged image');
    
    % Evaluate fitting
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
    
    % Save images
    imwrite(im2uint16(TrfImg), 'Reprojected.tif');
    imwrite(im2uint16(MergImg), 'Merged.tif');
    
    % Save evaluation of fitting
    fid = fopen('Fitting.txt', 'w');
    fprintf(fid, '%s\r\n', ['Number of pairs of points of interest: ' num2str(N)]);
    fprintf(fid, '%s\r\n', ['SSE, sum of squares due to error: ' num2str(SSE)]);
    fprintf(fid, '%s\r\n', ['MSE, mean squared error: ' num2str(MSE)]);
    fprintf(fid, '%s\r\n', ['RMSE, root mean squared error: ' num2str(RMSE)]);
    fprintf(fid, '%s\r\n', ['R_squared: ' num2str(R_squared)]);
    fprintf(fid, ['Magnitudes of residuals: ' num2str(ResidualMag')]);
    fclose(fid);