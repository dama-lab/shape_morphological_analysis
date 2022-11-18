function zscoreMatrix = zscoreReferenced(volMatrix, groupLabel, refGroup)
%% calculate the zscore of the matrix, with regards to the refGroup 
%  volMatix: raw matrix
%       row: subjects
%    column: structure/label/feature

% Calculating the zscore of the entire matrix

%% Determine reference group
refId = cellfun(@(x) strcmp(x,refGroup), groupLabel);
refGrpMatrix = volMatrix(refId,:);

%% Measure the mean/var of the reference group
refGrpMean = mean(refGrpMatrix,1);
refGrpVar = std(refGrpMatrix,0,1);

%% Calculate the zscore of volMatrix
zscoreMatrix = (volMatrix - refGrpMean)./refGrpVar;

end

