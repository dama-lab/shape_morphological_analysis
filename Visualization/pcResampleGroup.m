function [pcGroup, vertexValueOnMeanGroup, distMapGroup, dispMapGroup, ... 
    directDistMapGroup, vertexMean, pcMean] = pcResampleGroup(surfGroup, tformGroup, ...
    vertexValueDir, bgLabel, refTarget, groupMean)
    %% Groupwise Point-Cloud Registration
    %    surfGroup: containers.Map
    %    [Option 1]refSurf: reference surface (Recommented: from groupwise registration average)
    %    [Option 2]refTarget: target as reference surface

    %% by default, set background label to zero
    if ~exist('bgLabel','var') || isempty(bgLabel); bgLabel = 0; end
    if ~exist('refTarget','var') || isempty(refTarget); refTarget = 'groupMean'; end

    targetList = tformGroup.keys; % without 'groupMean'

    %% Get group mean surface vertex
    surfMean = surfGroup(refTarget); 
    
    %% Registration to the reference images
    pcGroup = containers.Map;
    % projected vertex value from each target onto group mean surface
    vertexValueOnMeanGroup = containers.Map; 
    % distance map for each target
    distMapGroup = containers.Map;
    % Directional Distance map (containing +/- directions)
    directDistMapGroup = containers.Map;
    % displacement map for each target;
    dispMapGroup = containers.Map;

    vertexRigidAll = [];
    for targetId = 1:length(targetList)
        target = targetList{targetId};
        fprintf('Resampling %d/%d: %s ...',targetId,length(targetList),target);
        
        %% skip reference target
        if ~isKey(tformGroup,target); continue; end

        %% get the surface vertex of the current target
        surfTarget = surfGroup(target);

        %% read/convert vertice to mask
        verticeNii = fullfile(vertexValueDir, target);
        vertice = single(niftiread(verticeNii));
        mask = (vertice ~= bgLabel);

        %% Find all point index of the mask
        maskIdx = find(mask);
        %% convert linear index into subscription
        [maskH,maskW,maskD] = ind2sub(size(mask),maskIdx);
        maskSub = [maskH, maskW, maskD];
        
        fprintf(' find vertex intensity from raw data ...');
        % for each vertex node in surfSub, get the index in the maskSub represending the
        % subscription of the nearest point on vol/maskSub (i.e. on the surface)
        [volSurfSub,~] = knnsearch(maskSub,surfTarget);


        %% get the vertex intensity value of the surface vertexs (closest points)
        % create point cloud only using the surface vertex with intensity from vol
        % get the subscription
        cmapSurfSub = maskSub(volSurfSub,:);
        % Convert subscription into idx for colormap
        cmapSurfIdx = sub2ind(size(vertice),cmapSurfSub(:,1),cmapSurfSub(:,2),cmapSurfSub(:,3));
        % get the vol intensity for each cmapSurIdx
        cmapSurf = vertice(cmapSurfIdx);

        %% convert surf to point cloud
        pcTarget = pointCloud(surfTarget,'intensity',cmapSurf);

        %% get the pre-computed rigid registration transform to reference (generated from ICP)   
        tformRigid = tformGroup(target);
        %% apply the transform to the point cloud
        pcTargetResampled = pctransform(pcTarget, tformRigid);
        pcGroup(target) = pcTargetResampled;
        %% find the closest value on the group mean surface/pointcloud
        fprintf(' project vertex intensity to group mean ...\n');
        % 1) find the corresponding points on transformed/registered/resampled surfTarget
        % for each vortex on surfMean 
        % 2) Also, get the Displacement/distance map from KNN search
        [surfLocSub,distMapGroup(target)] = knnsearch(pcTargetResampled.Location, surfMean);
        vertexValueOnMeanGroup(target) = pcTargetResampled.Intensity(surfLocSub);
        %% location of the closesest point
        vertexRigid = pcTargetResampled.Location(surfLocSub,:);
        vertexRigidAll = cat(3,vertexRigidAll,vertexRigid);
        %% Displacement map
        dispMapGroup(target) = vertexRigid - surfMean;

        %% determind Direction of distance (outside/within groupMean's surface)
        direction = int8(interp3(double(groupMean),vertexRigid(:,1),vertexRigid(:,2),vertexRigid(:,3),'spline'))>0;
        direction = (1-2*direction);
        directDistMapGroup(target) = direction .* distMapGroup(target);
    end
    
    vertexMean = mean(vertexRigidAll,3);
    pcMean = pointCloud(vertexMean);
end

