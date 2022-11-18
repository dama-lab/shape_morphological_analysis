function cortical_thickness_individual()
    %% calculate cortical thickness for all inividual subjects
    
    %% Load opts
%     plot_opt = load_opts;
%     v2struct(plot_opt);
%     

    rootDir = '/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/';
    %% Get Target List
    targetListDir = fullfile(rootDir,'TargetList');
    targetListWtFile = fullfile(targetListDir,'WT.txt');
    targetListTgFile = fullfile(targetListDir,'TG.txt');

    targetListWT = readTextLines(targetListWtFile);
    targetListTG = readTextLines(targetListTgFile);
    
    targetLists = [targetListWT, targetListTG];
        
    %%
    tic;
    thicknessDir = fullfile(rootDir,'06_thickness/groupMeanBackpropagate');
    for targetId = 22 %1:length(targetLists(:))
        % 22: tc1_274436-ob_c need reslicing
        %% Get target
        target = targetLists{targetId};
        fprintf('%s: ', target);
        %% load sublayer segemntation (WM1+GM2-Sulci)
        segVolNii = sprintf('%s/%s',WmGmRemoveSulciDir,target);
        rawVolNii = fullfile(rawImgDir,target);
        segVol = niftiread(segVolNii);

        %% Calculate Laplacian thickness (Full Cortex, Ganular, Molecular)
        [thickness, LaplacianField, L0, L1] = LaplacianThickness(segVol,2,1);
        %%
        % L0 = CSFD = molecular layer thickness
        % L1 = WMD = granular layer thickness
        % volPeek(L1);colormap jet
        
        %% Get nifti header info
        volNiftiInfo = niftiinfo(rawVolNii);
        volNiftiInfo.Datatype = 'double';
        
        %% Save thicknesses
        %% Save Full cortical thickness
        thicknessFullDir = fullfile(thicknessDir,'thickness_map');
        if ~exist(thicknessFullDir,'dir'); mkdir(thicknessFullDir); end
        niftiwrite(thickness, thicknessFile,volNiftiInfo,'Compressed',true); 
        
        %% Save Molecular layer thickness ( = CSFD = L0)
        thicknessMoleDir = fullfile(thicknessDir,'thickness_CSFD');
        if ~exist(thicknessMoleDir,'dir'); mkdir(thicknessMoleDir); end
        thicknessFile = fullfile(thicknessMoleDir, target);
        niftiwrite(L0, thicknessFile,volNiftiInfo,'Compressed',true);
        
        %% Save Granular layer thickness
        thicknessGranDir = fullfile(thicknessDir,'thickness_WMD');
        if ~exist(thicknessGranDir,'dir'); mkdir(thicknessGranDir); end
        thicknessFile = fullfile(thicknessGranDir, target);
        niftiwrite(L1, thicknessFile,volNiftiInfo,'Compressed',true);
        
    end
    toc
end

