% ImagePadding.m

% This script creates a black background image and places the input image
% at a position specified by the user on the aforementioned background.

% Written in 2020 by Buyun Tian @ Tao Xu's lab, Institute of Biophysics,
% Beijing, China.

clear;
close all;
clc;
[filename, filepath] = uigetfile('*.*','Select image to be padded');
InputImg = imread([filepath filename]);
Parameters = inputdlg({'New image width:', 'New image height:', ...
    'Postiion of original image (upper left x)', ...
    'Postiion of original image (upper left y)'}, ...
    'Parameters', 1, {'2048', '4096', '0', '0'}, 'on');
NewImgWidth = str2double(Parameters{1});
NewImgHeight = str2double(Parameters{2});
UpperLeft_x = str2double(Parameters{3});
UpperLeft_y = str2double(Parameters{4});
InputImgSize = size(InputImg);
InputImgWidth = InputImgSize(2);
InputImgHeight = InputImgSize(1);
PadTop = UpperLeft_x - 1;
PadBottom = NewImgHeight - UpperLeft_x - InputImgHeight + 1;
PadLeft = UpperLeft_y - 1;
PadRight = NewImgWidth - UpperLeft_y - InputImgWidth + 1;
Padded = padarray(InputImg, [PadTop PadLeft], 'pre');
Padded = padarray(Padded, [PadBottom PadRight], 'post');
imshow(Padded, []);
imwrite(Padded, 'Padded.tif');
