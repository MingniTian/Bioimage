function [centre_x, centre_y] = findbeadcentre2d(img)
%FINDBEADCENTRE 
%   [centre_x, centre_y] = findbeadcentre2d(img) finds the central position
%   (centre_x, centre_y) of a spherical bead on the input image matrix img
%   by matching disc-shaped templates with varying diameters to the image 
%   and choosing the best result.

% Written in 2020 by Buyun Tian @ Tao Xu's lab, Institute of Biophysics,
% Beijing, China.

imgsize = size(img);
if numel(imgsize) == 3
    img = rgb2gray(img); % Convert non-greyscale image to greyscale
    imgsize = size(img);
end
img = imadjust(1 - img); % Make bead appear bright
dmin = max(3, floor(0.2*min(imgsize))); % Diameter lower bound
dmax = min(imgsize); % Diameter upper bound
corrmax = zeros(1,dmax); % Maximum correlation value for every diameter
peakpositions = zeros(dmax, 2); % Position of correlation peak for every diameter
for d = dmin:dmax % Varying diameters
    mask = gendiscmask(d); % Disc-shaped mask with diameter d
    masksize = size(mask);
    c = normxcorr2(mask, img); % Normalised cross-correlation
    c = normxcorr2post(c, masksize); % Values computed with zero-padded edges set to -1
    corrmax(d) = max(max(c));
    [ypeak, xpeak] = find(c == corrmax(d));
    xpeak = mean(xpeak);
    ypeak = mean(ypeak);
    peakpositions(d,1) = xpeak - (masksize(2) - 1)/2;
    peakpositions(d,2) = ypeak - (masksize(1) - 1)/2;
end
bestd = find(corrmax == max(max(corrmax))); % Best matched diameter
centre_x = peakpositions(bestd, 1);
centre_y = peakpositions(bestd, 2);
end
