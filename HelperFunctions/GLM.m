function [volStdResidMtx,volResidMtx] = GLM(volMtx, covar,groupLabel, refGroup)
%% calculate the GLM of the volMtx, return only the residual
%  by regressing out the effiect of covariate (currently only for TIV)
%% Input parameters
%  volMtx:    raw matrix
%       row:    subjects
%       column: vertex/structure/label/feature
%  covar: covariate to be regressed out e.g. TIV
%       row:  = no of subjects
%       colume: individual covariates

%% Extract control group only:
refId = cellfun(@(x) strcmp(x,refGroup), groupLabel);
grpMtxRef = volMtx(refId,:);
covarRef = covar(refId==1,:);

%% GLM model fitting (take about 3 minutes for 10000 points)
tic
volResidMtx = nan(size(volMtx));
for vertex = 1:size(grpMtxRef,2)
    if mod(vertex,1000) == 0; fprintf('%d ',vertex); end
    featureRaw = volMtx(:,vertex);
    featureRawRef = grpMtxRef(:,vertex);
    mdl = fitlm(covarRef,featureRawRef,'linear');
    volPred = predict(mdl,covar);
    residual = featureRaw - volPred;
    volResidMtx(:,vertex) = residual;
    if mod(vertex,10000) == 0
        fprintf('\n')
    end
end
toc

%% Create w-score (z-score of the residual matrix)
volStdResidMtx = zscoreReferenced(volResidMtx,groupLabel, refGroup);


%% %%%%% mvregress
% %% n: number of samples (observations), d: number of vertex (measurements)
% [n,d] = size(volMatrix);
% Xmat = [ones(n,1), covar];
% 
% %% Multivariate linear regression
% % Ref: https://www.mathworks.com/help/stats/mvregress.html#btkr6ds-3
% % Example: https://www.mathworks.com/help/stats/multivariate-general-linear-model.html
% 
% % design matrix
% Xcell = cell(1,n);
% % kro: kronecker tensor product
% for i = 1:n
%     Xcell{i} = kron(Xmat(i,:),eye(d));
% end
% 
% [beta,sigma,residualMatrix,V] = mvregress(Xcell,volMatrix);
% residualMatrix;


