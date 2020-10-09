function [centre_x, centre_y] = gausscentre2d(img)
%GAUSSCENTRE2D
%   [centre_x, centre_y] = gausscentre2d(img) finds the central position
%   (centre_x, centre_y) of a spot on the input image matrix img by 2D
%   Gaussian fitting.

%   Written in 2020 by Buyun Tian @ Tao Xu's lab, Institute of Biophysics,
%   Beijing, China.
imgsize = size(img);
if numel(imgsize) == 3
    img = rgb2gray(img); % Convert non-greyscale image to greyscale
    imgsize = size(img);
end
BG_est = (sum(img(1:imgsize(1),1)) + sum(img(1:imgsize(1),imgsize(2))) ...
         + sum(img(1,2:(imgsize(2) - 1))) ...
         + sum(img(imgsize(1),2:(imgsize(2) - 1))))...
         /((imgsize(1)+imgsize(2))*2 - 4); % Estimate background intensity
x = ones(imgsize(1), 1)*(1:imgsize(2));
y = (1:imgsize(1))'*ones(1, imgsize(2));
sumintensity = sum(sum(img - BG_est));
centre_x_est = sum(sum((img - BG_est).*x))/sumintensity; % Estimate central x
centre_y_est = sum(sum((img - BG_est).*y))/sumintensity; % Estimate central y
clear x y;
I_est = max(max(img)) - BG_est; % Estimate maximum intensity
sigma_est = sqrt(sumintensity/(2*pi*I_est)); % Estimate gaussian sigma
[x, y, value] = prepareSurfaceData(1:imgsize(2), 1:imgsize(1), img);
gauss2d = fittype('BG + I*exp(-((x - x0)^2 + (y - y0)^2)/(2*sigma^2))',...
                  'coefficients', {'BG', 'I', 'x0', 'y0', 'sigma'},...
                  'dependent', {'value'}, 'independent', {'x', 'y'});
foptions = fitoptions(gauss2d);
foptions.StartPoint = [BG_est I_est centre_x_est centre_y_est sigma_est];
foptions.Lower = [0 0 0 0 0];
foptions.Upper = [max(max(img)) 1.5*max(max(img)) ...
                  imgsize(2)+1 imgsize(1)+1 min(imgsize)/2];
f = fit([x, y], value, gauss2d, foptions);
centre_x = f.x0;
centre_y = f.y0;
end