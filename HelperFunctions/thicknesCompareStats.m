function [pc, tstats, volStdResidMtx, vertexValueOnMeanGroup] = thicknesCompareStats(targetLists, TIVcsv, groupMean, ...
    surfGroup, tformGroup, vertexValueDir, bgLabel, resultMat)
    
    %% Resample surface to group mean
    [pcGroup,vertexValueOnMeanGroup,distMapGroup,dispMapGroup,directDistMapGroup, vertexMean, pcMean] ...
        = pcResampleGroup(surfGroup, tformGroup, vertexValueDir, bgLabel, 'groupMean', groupMean);
    
    %% compute groupwise statistics
    % get the TIV
    TIVtable = readtable(TIVcsv);
    TIV = TIVtable.('TIV');
    groupLabel = TIVtable.('population');
    refGroup = "wt";
    % smooth the resampled vertex value   
    surfMatrix = group2Matrix(targetLists,vertexValueOnMeanGroup);
    surfMatrixSmooth = smoothdata(surfMatrix,2,'gaussian',10);
    
    %% GLM-based TIV
    [volStdResidMtx,volResidMtx] = GLM(surfMatrixSmooth, TIV, groupLabel, refGroup);
    
    %% unpaired t-test with FDR=0.05
    % Get group index
    volNormalize = volStdResidMtx;
    wt_idx = (groupLabel==refGroup);
    tg_idx = ~(wt_idx);
    for vertex = 1:size(volNormalize,2)
        % get group data
        volTW = volNormalize(wt_idx,vertex);
        volTG = volNormalize(tg_idx,vertex);
        [tstats.h(vertex),tstats.p(vertex),tstats.ci(vertex,:),tstats.stats(vertex)] = ttest2(volTW,volTG);
    end
    %% Multiple Comparison Correction with False Discovery Rate (FDR) =0.05
    fdr_q = 0.1;
    [tstats.adj_h, tstats.crit_p, tstats.adj_ci_cvrg, tstats.adj_p]=fdr_bh(tstats.p,fdr_q,'pdep','yes');
    
    %% create t-statistical map
    refSurf = surfGroup('groupMean');
    % adjust the displayed p-value
    adjPshow = smoothdata(tstats.adj_p'); % smooth
    adjPshow(tstats.adj_h==0) = 1;
    adjPshow = 1 - adjPshow;
    % show on PtCloud
    pc.tstats = pointCloud(refSurf,'Intensity',adjPshow);
    
    %% create mean cortical distance for two groups
    volNormalize = volStdResidMtx;
    % Get group index
    wt_idx = (groupLabel==refGroup);
    tg_idx = ~(wt_idx);
    % get group volume mean
    volWT = mean(volNormalize(wt_idx,:),1);
    volTG = mean(volNormalize(tg_idx,:),1);
    
    pc.diff = pointCloud(refSurf,'Intensity',(volTG' - volWT'));    
    
    
    %% Save variables
    save(resultMat,'pcGroup','vertexValueOnMeanGroup','distMapGroup','directDistMapGroup','pc','tstats');
end


