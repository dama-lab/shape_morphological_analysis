function [slm, stat] = GLMSurfStats(surfThickMat, surfaceSurfstat, covars, contrast_var)
    %% GLM using Surfstats
    
    %% Define formular (covar1 = TIV; covar2 = group)
    Model = 1;
    covarNames = fieldnames(covars);
    for covarNo = 1:numel(covarNames)
        covarName = covarNames{covarNo};
        covar = covars.(covarName);
        Model = Model + term(covar,covarName);
    end
    % figure; image(Model)
    
    %% setup linear model
    slm = SurfStatLinMod(permute(surfThickMat,[2,1,3]), Model, surfaceSurfstat);
    
    %% group contrast
    if ~exist('contrast_var','var') || isempty(contrast_var)
        contrat_var = "group"; % for backward compatible
    end
    Group = term(covars.(contrast_var));
    contrast = Group.wt - Group.tg;
    
    slm = SurfStatT(slm, contrast);
    
    %% [Visualize] T-statistics
    % figure;SurfStatView(slm.t, surfaceStat, 'T')
    
    %% [Test] (optional) To find the threshold for P=0.05, corrected, the resels are:
    stat.resels = SurfStatResels(slm);
    % stat_threshold(resels, length(slm.t), 1, slm.df);
    
    %% [Stats] However the best way is to view the P-values for each vertex.
    % pval.P contains P-values for peaks, and pval.C contains P-values for clusters.
    [stat.pval, stat.peak, stat.clus] = SurfStatP(slm);
    
    %% Calculate Q-value
    stat.qval = SurfStatQ(slm);
    
    
    %% [Visualize] Model
    % This special structure of pval is recognised by SurfStatView which draws the figure in a special way:
    % figure; SurfStatView(pval, surfaceStat, 'WT v.s. TG removing ICV')
    