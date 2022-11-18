function sub = mask2sub(vol, bgLabel)
    %% find the subscription index of all the foreground voxels in mask
    
    %% convert vol to mask ([0,1])
    if ~exist('bgLabel','var') || isempty(bgLabel)
        bgLabel = 0; 
    end
    vol = (vol ~= bgLabel);
    
    %% Find all point index of the mask
    maskIdx = find(vol);

    %% convert linear index into subscription
    [maskH,maskW,maskD] = ind2sub(size(vol),maskIdx);
    sub = [maskH, maskW, maskD];
    
    %% Find all point index of the mask
    maskIdx = find(vol);
    
    %% convert linear index into subscription
    [maskH, maskW, maskD] = ind2sub(size(vol), maskIdx);
    sub = [maskH, maskW, maskD];