function groupMtx = group2Matrix(targetLists, vertexValues)
%% construct vertex-sise matrix to prepare for GLM
groupMtx = [];
for targetId = 1:length(targetLists(:))
    target = targetLists(targetId);
    %% skip last empty target
    if ~strlength(target); break; end
    %% read vertex-wise value
    groupMtx = [groupMtx; vertexValues(target)'];
end
