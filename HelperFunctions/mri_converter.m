function opt = mri_converter(opt)
    %% %%%%%%%%%%
    % MRI volume/surface/mesh/point-cloud conversion
    % Definition:
    %      layer: skeletonized layer (e.g. Purkinje layer)
    %      pc: point cloud
    %% %%%%%%%%%%

    if strcmp(opt.convert_type, 'layer2pc') % no need anymore
        opt.pc = layer2pc(opt.intensity,opt.layer);
    elseif strcmp(opt.convert_type, 'vol2pc')
        % generate mask from vol if does not exist
        if ~isfield(opt, 'mask'); opt.mask = opt.vol; end
        if ~isfield(opt, 'bgLabel'); opt.bgLabel = 0; end
        if isfield(opt, 'iso2meshOpt')
            opt.pc = vol2pc(opt.vol, opt.mask, opt.bgLabel, opt.iso2meshOpt);
        else
            opt.pc = vol2pc(opt.vol, opt.mask, opt.bgLabel);
        end
    elseif strcmp(opt.convert_type,'vol2surf')
        [opt.surf,opt.elem,opt.regions, opt.holes] = vol2surf(opt.mask, opt.iso2meshOpt, opt.bgLabel);
    end
end % end of main function

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individual functions
function pc = layer2pc(intensity,layer)
    %% convert 3d thin layer to 3d point cloud, with intensity colormap
    % specified in intensity
    
    % point cloud index
    idx = find(layer);
    % convert to x-y-z- coordinator
    [x, y, z] = ind2sub(size(layer),idx);
    % color on the pc
    pcColor = intensity(idx);
    % create point cloud
    pc = pointCloud([x,y,z],'Intensity',pcColor);
end

function pc = vol2pc(vol,mask, bgLabel, iso2meshOpt)
    %% convert 3d vol into point cloud, with intensity  value assigned as the surface value
    % vol: 3D volumetric matrix contain intensity colormap info (i.e. struct-label, volume, thickness)
    % [Optional]: mask used to create point cloud (default = vol) has to be the same size with vol
    % [Optional]: iso2meshFlag (default = true ie. logical 1)
    %% if mask not provided, use the non-zero value of the vol
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
    maskIdx = find(mask);
    
    %% convert linear index into subscription
    [maskH,maskW,maskD] = ind2sub(size(mask),maskIdx);
    maskSub = [maskH, maskW, maskD];
    
    if iso2meshOpt.iso2meshFlag==true
        %% Recommended option: construct point cloud using only points from mask
        % the iso2mesh take some extra time to generate, but only use
        % surface points to create point cloud, so memory/cpu efficient

        %% define default iso2mesh parameters
        if ~isfield(iso2meshOpt,'maxnode'); iso2meshOpt.maxnode = 10000; end
        if ~isfield(iso2meshOpt,'radbound'); iso2meshOpt.radbound = 2; end
        if ~isfield(iso2meshOpt,'dofix'); iso2meshOpt.dofix=1; end
        if ~isfield(iso2meshOpt,'isovalues'); iso2meshOpt.isovalues=0.5; end
        if ~isfield(iso2meshOpt,'method'); iso2meshOpt.method = 'cgalsurf'; end % 'cgalmesh' takes 3x more time (20 sec vs  60 sec)
        % otherwise, simply use: opt = 0.5;
        
        %% find surfSub using v2s
        [surfSub,elem,regions,holes] = v2s(uint8(mask>0),iso2meshOpt.isovalues, iso2meshOpt, iso2meshOpt.method);
        
        % for each vertex node in surfSub, get the index in the maskSub represending the
        % subscription of the nearest point on vol/maskSub (i.e. on the surface)
        [volSurfSub,~] = knnsearch(maskSub,surfSub);
        
        % get the colormap intensity value of the surface vertexs (closest points)
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


function [surfSub,elem,regions,holes] = vol2surf(mask, iso2meshOpt, bgLabel)
    %% convert 3d vol into point cloud, with intensity  value assigned as the surface value
    % vol: contain intensity colormap info (i.e. struct-label, volume, thickness)
    % [Optional]: mask used to create point cloud (default = vol) has to be the same size with vol
    % [Optional]: iso2meshFlag (default = true ie. logical 1)
    %% if mask not provided, use the non-zero value of the vol
    if ~exist('iso2meshOpt','var')
        iso2meshOpt.iso2meshFlag = true;
    end
    
    %% convert mask to [0,1] (in case of 0,255)
    mask = (mask ~= bgLabel);
    
    %% Find all point index of the mask
    maskIdx = find(mask);
    
    %% convert linear index into subscription
    [maskH,maskW,maskD] = ind2sub(size(mask),maskIdx);
    maskSub = [maskH, maskW, maskD];

    %% construct surface point coordination using only points from mask surface
    % the iso2mesh take some extra time to generate, but only use
    % surface points to create point cloud, so memory/cpu efficient

    %% define default iso2mesh parameters
    if ~isfield(iso2meshOpt,'maxnode'); iso2meshOpt.maxnode = 10000; end % 40000
    if ~isfield(iso2meshOpt,'radbound'); iso2meshOpt.radbound = 2; end
    if ~isfield(iso2meshOpt,'dofix'); iso2meshOpt.dofix=1; end
    if ~isfield(iso2meshOpt,'isovalues'); iso2meshOpt.isovalues=0.5; end
    if ~isfield(iso2meshOpt,'method'); iso2meshOpt.method = 'cgalsurf'; end % 'cgalmesh' takes 3x more time (20 sec vs  60 sec)
    % otherwise, simply use: opt = 0.5;

    %% find surfSub using v2s
    [surfSub,elem,regions,holes] = v2s(uint8(mask>0),iso2meshOpt.isovalues, iso2meshOpt, iso2meshOpt.method);

end
