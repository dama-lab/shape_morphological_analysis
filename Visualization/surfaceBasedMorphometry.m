function residualMatrix = surfaceBasedMorphometry(groupMtx, TIVtable,refGroup)
%% based on the displacement field (or distance between each surface)
% To be implemented later
%   Detailed explanation goes here

%% Calculate GLM to remove the effect of TIV and only keep the residual
TIV = TIVtable.('TIV');
groupLabel = TIVtable.('population');
% Multivariate linear regression Ref: https://www.mathworks.com/help/stats/mvregress.html#btkr6ds-3
% Example: https://www.mathworks.com/help/stats/multivariate-general-linear-model.html
residualMatrix = GLM(groupMtx, TIV,groupLabel, refGroup);

%% compute the Wscore ([Optional] for better visual representationonly, no effect on the statisticals)
zscoreReferenced


end

