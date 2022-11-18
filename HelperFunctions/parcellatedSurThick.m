function meanThicknesses = parcellatedSurThick(vertexThickness, vertexParcellation, labels)
    %% parcellatedSurThick(vertexThickness, vertexParcellation, labels)
    % get mean surface thickness for each parcellated label from vertex information
    
    %% Automatically determine label numbers is not already specified
    if ~exist('labels','var') || isempty(labels)
        labels = unique(vertexParcellation);
    end
    
    %% calculate vertex mean thickness
    meanThicknesses = nan(size(labels));
    labelNo = 1;
    for label = labels
        meanThicknesses(labelNo) = mean(vertexThickness(vertexParcellation == label));
        labelNo = labelNo + 1;
    end