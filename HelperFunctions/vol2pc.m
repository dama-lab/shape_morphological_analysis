function pc = vol2pc(vol,mask, bgLabel, iso2meshOpt)
    %% convert 3d vol into point cloud, with intensity  value assigned as the surface value
    % vol: 3D volumetric matrix contain intensity colormap info (i.e. struct-label, volume, thickness)
    % [Optional]: mask used to create point cloud (default = vol) has to be the same size with vol
    % [Optional]: iso2meshFlag (default = true ie. logical 1)
    %% if mask not provided, use the non-zero value of the vol
    if ~exist('bgLabel','var') || isempty(bgLabel)
        bgLabel = 0; 
    end
    if ~exist('mask','var'); mask = vol; end
    if ~exist('iso2meshOpt','var')
        iso2meshOpt.iso2meshFlag = true;
    end
    
    %% convert mask to [0,1] (in case of 0,255)
    mask = (mask ~= bgLabel);
    
    %% check if vol and mask have same size
    if sum(~(size(vol) == size(vol)))
        errir('vol2pc: vol and mask size mismatch');
        return
    end
    
    %% Find all point index of the mask
    % maskSub = mask2sub(mask, bgLabel);
    
    maskIdx = find(mask);
    
    %% convert linear index into subscription
    [maskH,maskW,maskD] = ind2sub(size(mask),maskIdx);
    maskSub = [maskH, maskW, maskD];
    
    if iso2meshOpt.iso2meshFlag==true
        %% Recommended option: construct point cloud using only points from mask
        % the iso2mesh take some extra time to generate, but only use
        % surface points to create point cloud, so memory/cpu efficient

        %% define default iso2mesh parameters
        % iso2meshOot = iso2meshDefaultOpt()
        if ~isfield(iso2meshOpt,'maxnode'); iso2meshOpt.maxnode = 40000; end
        if ~isfield(iso2meshOpt,'radbound'); iso2meshOpt.radbound = 2; end
        if ~isfield(iso2meshOpt,'dofix'); iso2meshOpt.dofix=1; end
        if ~isfield(iso2meshOpt,'isovalues'); iso2meshOpt.isovalues=0.5; end
        if ~isfield(iso2meshOpt,'method'); iso2meshOpt.method = 'cgalsurf'; end % 'cgalmesh' takes 3x more time (20 sec vs  60 sec)
        % otherwise, simply use: opt = 0.5;
        
        %% find surfSub using v2s
        [surfSub,elem,regions,holes] = v2s(uint8(mask>0),iso2meshOpt.isovalues, iso2meshOpt, iso2meshOpt.method);
        
        %% for each vertex node in surfSub, get the index in the maskSub represending the
        % subscription of the nearest point on vol/maskSub (i.e. on the surface)
        [volSurfSub,~] = knnsearch(maskSub,surfSub);
        
        %% get the colormap intensity value of the surface vertexs (closest points)
        % get the subscription
        cmapSurfSub = maskSub(volSurfSub,:);
        % Convert subscription into idx for colormap
        cmapSurfIdx = sub2ind(size(vol),cmapSurfSub(:,1),cmapSurfSub(:,2),cmapSurfSub(:,3));
        % get the vol intensity for each cmapSurIdx
        cmapSurf = vol(cmapSurfIdx);
        % create point cloud only using the surface vertex with intensity from vol
        pc = pointCloud(surfSub,'intensity',cmapSurf);
    else
        %% Alternative option: construct point cloud using all points inside mask
        % memory/cpu inefficient, not recommented
        % get the intensity value of all the points inside the vol(linear indexed)
        volIntensity = vol(maskIdx);
        % create point cloud using all vol points inside mask
        pc = pointCloud(maskSub,'Intensity',volIntensity);
    end
end