function [surfDownsampled,pcDownsampled] = downsampleSurface(surface,downsampleRate)
%% Downsample surf by internally converting the into point cloud
    % Convert surf to pc
    pc = pointCloud(surface);
    % downsample point cloud
    pcDownsampled = downsamplePc(pc,downsampleRate);
    % extract downsampled surf from point cloud
    surfDownsampled = pcDownsampled.Location;
end

