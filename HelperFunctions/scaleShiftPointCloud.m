function scaleShiftPointCloud
    %% Scale and shift a point cloud
    % Ref: https://www.mathworks.com/matlabcentral/answers/71458-3d-registration-for-2-clouds-of-points
    
    % Do this for all three axises
    range = max(x)-min(x);
    normalizedX = desiredRange * (x-mean(x)) / range + desiredMean;
    
end