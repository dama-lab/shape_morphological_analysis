function volSize =  extractVol(segVol,labels,voxeldim)
    %% Extract volumes (voxel number) of labels in the segVol
    
    %% Automatically determine label numbers is not already specified
    if ~exist('labels','var') || isempty(labels)
        labels = unique(segVol)';
    end
    %% Automatically set voxel size = 1 is not already specified
    if ~exist('voxeldim','var') || isempty(labels)
        voxeldim = [1,1,1];
    end
    
    %% count voxel numbers
    volSize = nan(size(labels));
    labelNo = 1;
    for label = labels
        % First: calculate voxelSum
        volume = sum(segVol(:) == label);
        % times pixel/voxel dimension
        for dim = 1:length(voxeldim)
            volume = volume * voxeldim(dim);
        end
        volSize(labelNo) = volume;
        labelNo = labelNo + 1;
    end
    
    
    