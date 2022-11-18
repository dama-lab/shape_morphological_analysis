function ExtractBiggestBlob()
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;
% Read in a standard MATLAB gray scale demo image.
folder = fullfile(matlabroot, '\toolbox\images\imdemos');
baseFileName = 'coins.png';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
  % File doesn't exist -- didn't find it there.  Check the search path for it.
  fullFileName = baseFileName; % No path this time.
  if ~exist(fullFileName, 'file')
    % Still didn't find it.  Alert user.
    errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
    uiwait(warndlg(errorMessage));
    return;
  end
end
grayImage = imread(fullFileName);
% Get the dimensions of the image.
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(grayImage);
% Display the original gray scale image.
subplot(2, 2, 1);
imshow(grayImage, []);
title('Original Grayscale Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Give a name to the title bar.
set(gcf,'name','Demo by ImageAnalyst','numbertitle','off')
% Let's compute and display the histogram.
[pixelCount, grayLevels] = imhist(grayImage);
subplot(2, 2, 2);
bar(pixelCount);
grid on;
title('Histogram of original image', 'FontSize', fontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.
% Threshold the image to binarize it.
binaryImage = grayImage > 100;
% Fill holes
binaryImage = imfill(binaryImage, 'holes');
% Display the image.
subplot(2, 2, 3);
imshow(binaryImage, []);
title('Binary Image', 'FontSize', fontSize);
% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
[labeledImage, numberOfBlobs] = bwlabel(binaryImage);
blobMeasurements = regionprops(labeledImage, 'area', 'Centroid');
% Get all the areas
allAreas = [blobMeasurements.Area] % No semicolon so it will print to the command window.
menuOptions{1} = '0'; % Add option to extract no blobs.
% Display areas on image
for k = 1 : numberOfBlobs           % Loop through all blobs.
  thisCentroid = [blobMeasurements(k).Centroid(1), blobMeasurements(k).Centroid(2)];
  message = sprintf('Area = %d', allAreas(k));
  text(thisCentroid(1), thisCentroid(2), message, 'Color', 'r');
  menuOptions{k+1} = sprintf('%d', k);
end
% Ask user how many blobs to extract.
numberToExtract = menu('How many do you want to extract', menuOptions) - 1;
% Ask user if they want the smallest or largest blobs.
promptMessage = sprintf('Do you want the %d largest, or %d smallest, blobs?',...
  numberToExtract, numberToExtract);
titleBarCaption = 'Largest or Smallest?';
sizeOption = questdlg(promptMessage, titleBarCaption, 'Largest', 'Smallest', 'Cancel', 'Largest');
if strcmpi(sizeOption, 'Cancel')
  return;
elseif strcmpi(sizeOption, 'Smallest')
  % If they want the smallest, make the number negative.
  numberToExtract = -numberToExtract;
end