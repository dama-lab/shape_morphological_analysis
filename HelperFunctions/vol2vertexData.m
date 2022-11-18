function vertexIntensity = vol2vertexData(vertexSub, vol, mask, bgLabel, interp)
    %% find the intensity(colormap) of vertex points on volsurfSub
    % Input:
    %   vertexSub:  vertex subscription of the enquiring surface
    %   vol:        volume containing the intensity value of the surface
    %   mask:       (optional) mask to constrain the region/realm to
    %                          compute the interpolated vertex intensityvertexSub
    %   bgLabel: 
    %   interp: number of closest points to interpolate.
    %       if multiple closes points are included, calculate the weighted
    %       mean using the average of n closest points on mask
    %       Recommended value:
    %       0): directly interpolate from the entire 3D volumetric field
    %           with "interp3",(useful for surface displacement calculation)
    %       1) 1 closest point: = nearest neighbout interpolation
    %       4) 4 closes points: ~ project + in-plane 2D linear interpolation outside mask region
    %       6) 6 closes points: ~ balanced inside/outside mask region
    %       8) 8 closes points: ~ 3D linear interpolation inside mask region
    %   interpMethod: [Optional] surface point interpolation methods
    %       
    %            

    
    %% convert mask to [0,1] (in case of 0,255)
    if ~exist('bgLabel','var') || isempty(bgLabel)
        bgLabel = 0; 
    end
    if ~exist('mask','var') || isempty(mask)
        mask = (vol > bgLabel); 
    end
    % Using linear interpolation by default
    if ~exist('interp','var') || isempty(interp)
        interp = 1; 
    end
    
    %% Direct interpolation (doesn't work for thickness, maybe good for displacement map)
    if interp == 0
        vertexIntensity = interp3(vol,vertexSub(:,1),vertexSub(:,2),vertexSub(:,3),'makima');
        return;
    end
    %% find all the subscription index from mask
    maskSub = mask2sub(vol, bgLabel);
    
    %% for each vertex node in surfSub, get the index in the maskSub represending the
    % subscription of the nearest point on vol/maskSub (i.e. on the surface)
    [knnMaskSub,vertexDist] = knnsearch(maskSub, vertexSub, 'K', interp);
    %% Calculate the interpolation weights based on the vertedDist
    interpWeights = vertexDist ./ sum(vertexDist,2);
    %% Get the intensity subscription of k-nearest points
    volIntensity = [];
    for pointId = 1:size(knnMaskSub,2)
        %% Get the intensity subscription of k-nearest points
        IntSub = maskSub(knnMaskSub(:,pointId),:);
        % intensitySub(:,pointId,:) = IntSub;
        %% Convert sub-to-index
        IntInd = sub2ind(size(vol),IntSub(:,1),IntSub(:,2),IntSub(:,3));
        %% Get the intensity value of each of the k-nearest points
        volIntensity(:,pointId) = vol(IntInd);
    end
    %% Linear interpolate the vertex intensity colormap data
    vertexIntensity = dot(volIntensity,interpWeights,2);


    
    
    
    
    
    