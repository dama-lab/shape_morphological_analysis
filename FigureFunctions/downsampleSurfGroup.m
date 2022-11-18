function [pcDownsampleGroup,surfDownsampleGroup] = downsampleSurfGroup(targetLists, downsampleRate,vertexValueOnMeanGroup,refSurf)
%% downsample/resample (by converting to ptCloud first)
% Ref: https://www.mathworks.com/help/vision/ref/pcdownsample.html#bupqqn1-1-gridStep
%
% resampleRate: float number. <1: downsample; >1 upsample; =1: no change
% vertexValueOnMeanGroup: i.e. thickness(used to generate ptCloud)
% refSurf: surface coordinate of group mean [1000x3]
    pcDownsampleGroup = containers.Map;
    surfDownsampleGroup = containers.Map;
    %%
    for targetId =1:length(targetLists(:))
        target = targetLists(targetId);
        
        %% skip empty target
        if ~strlength(target); continue; end
        
        %% extract the vertex value for the current target
        vertexValueTarget = vertexValueOnMeanGroup(target);
        
        %% convert surface+intensity to ptCloud
        pcTarget = pointCloud(refSurf,'Intensity',vertexValueTarget);

        %% resample pointcloud into evenly distributed grid
        % pcTargetResample = pcdownsample(pcTarget,'gridAverage',0.005);
        
        %% downsample pointcloud
        stepSize = round(1/downsampleRate);
        indices = 1:stepSize:pcTarget.Count;
        pcDownsampleGroup(target) = select(pcTarget,indices);
        surfDownsampleGroup(target) = pcDownsampleGroup(target).Intensity;
    end

end