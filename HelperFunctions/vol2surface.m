function surface = vol2surface(vol,mask, bgLabel, iso2meshOpt, mask2subFlag,interp)
    %% convert 3d vol into point cloud, with intensity  value assigned as the surface value
    % mask used to create point cloud (default = vol) has to be the same size with vol
    % [Optionsl]: bgLabel - Background label number. Default = 0 (normally 0 or NaN)
    % [Optional]: iso2meshFlag (default = true ie. logical 1)
    % [Optional]: interp: number of closest points to interpolate
    %% if mask not provided, use the non-zero value of the vol
    % 
    
    if ~exist('bgLabel','var') || isempty(bgLabel)
        bgLabel = 0; 
    end
    if ~exist('mask','var') || isempty(mask)
        mask = (vol > bgLabel); 
    end
    if ~exist('iso2meshOpt','var') || isempty(iso2meshOpt)
        iso2meshOpt.iso2meshFlag = true;
    end
    if ~exist('mask2subFlag','var'); mask2subFlag=0; end
    if ~exist('interp','var') || isempty(mask)
        interp = 1; 
    end
    %% %%%%%%%%%%% Alternative way: getting pointcloud first, then convert pc to surface 
    % (Better choice for very thin mask, although results are inferior compared to the default option using iso2mesh)
    if mask2subFlag == 1
        %% convert mask to [0,1] (in case of 0,255)
        mask = uint8(mask ~= bgLabel);
       
        
        %% Find all point index of the mask
        maskIdx = find(mask);

        %% convert linear index into subscription
        % Can upsample first, and use 0.5 threshold
        [maskH,maskW,maskD] = ind2sub(size(mask),maskIdx);
        clear surface
        surface.vertex = [maskH, maskW, maskD];
        % construct triedge from vertex points, yet to implement
        pt2surf(surface.vertex);
    else
        %% construct surface point coordination using only points from mask surface
        % the iso2mesh take some extra time to generate, but only use
        % surface points to create point cloud, so memory/cpu efficient

        %% define default iso2mesh parameters
        % iso2meshOot = iso2meshDefaultOpt()
        if ~isfield(iso2meshOpt,'maxnode'); iso2meshOpt.maxnode = 20000; end % 40000
        if ~isfield(iso2meshOpt,'radbound'); iso2meshOpt.radbound = 2; end
        if ~isfield(iso2meshOpt,'dofix'); iso2meshOpt.dofix=1; end
        if ~isfield(iso2meshOpt,'isovalues'); iso2meshOpt.isovalues=0.5; end
        if ~isfield(iso2meshOpt,'method'); iso2meshOpt.method = 'cgalsurf'; end % 'cgalmesh' takes 3x more time (20 sec vs  60 sec)
        % otherwise, simply use: opt = 0.5;

        %% find surfSub using v2s
        [surface.vertex, surface.triedge, surface.regions, surface.holes] = ...
            v2s(mask,iso2meshOpt.isovalues, iso2meshOpt, iso2meshOpt.method);

        %% fix mesh generated by iso2mesh by removing isolated nodes (duplicated non-manifold vertices) (Ref: AnalyzeToParaviewVtk.m)
        [surface.vertex, surface.triedge] = meshcheckrepair(surface.vertex, surface.triedge);
        
        %% Get the surface intensity
        surface.intensity = vol2vertexData(surface.vertex, vol,[],[],interp);
    end
