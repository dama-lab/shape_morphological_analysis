function cortical_sublayer_thickness(plot_opt)
    %% load all variables
    v2struct(plot_opt);
    if ~exist(thicknessGranularDir,'dir'); mkdir(thicknessGranularDir); end 
    if ~exist(thicknessMoleDir,'dir'); mkdir(thicknessMoleDir); end 

    %% 
    targetNo = length(targetLists(:));
    
    for targetId = 1%:targetNo
        target = targetLists(targetId);
        % target = 'tc1_275322-ob_c.nii.gz';
        % skip empty target (the last one)
        if strcmp(target,""); continue; end
        fprintf('\n%s: \n',target);
        
        equiVolumeNii = fullfile(thickness)
        
        %% %%%%%%%%%%% Alternative approach %%%%%%%%%%%%%%%
        % calculate all thickness on granular layer
        % Doesn't seems to work well with surface plot
        
        %% define file locations
        granularWmNii = fullfile(granularWmNiiDir,target);
        thicknessFullNii = fullfile(thicknessCorticalFullDir,target);
        thicknessGranNii = fullfile(thicknessGranularDir,target);
        thicknessMoleNii = fullfile(thicknessMoleDir, target);
        
        %% Read Granular (1) + cortical_sublayer_thickness(plot_opt)WM (2)
        granularWm = niftiread(granularWmNii);
        % volPeek(granularWm);colormap jet
        
        %% Calculate Laplacian thickness
        % thickness = LaplacianThickness(granularWm, 1, 2);
        fprintf(' Granular LaplacianField ... \n');
        LalpacianField = LaplaceField(granularWm, 1, 2, 1);
        %%
        fprintf(' Granular Thickness ... \n');
        [thicknessGran, L0, L1, epsilon_array] = laplace_thickness(LalpacianField);

        %% Save granular nifti file
        niftiwrite(thicknessGran, thicknessGranNii, ...
            niftiinfo(thicknessFullNii), 'Compressed',true);
%         
        %% calculate molecular thickness
        fprintf(' Molecular Thickness ... ');
        %% Option 1: simply substract from full thickness
        thicknessGran = niftiread(thicknessGranNii);
        thicknessFull = niftiread(thicknessFullNii);
        thicknessMole = thicknessFull - thicknessGran;
        thicknessMole(thicknessGran==0)=0; 
        %volPeek(thicknessMole); colormap jet
        
        %% Option 2: calculate within molecular region
%         gmRemoveSulciNii = fullfile(gmRemoveSulciDir, target, 'GM_with_sulci.nii.gz');
%         gmRemoveSulci = niftiread(gmRemoveSulciNii);
%         % volPeek(gmRemoveSulci);
%         %% extract molecular mask
%         moleGran = gmRemoveSulci;
%         moleGran(granularWm==1)=2;
%         %volPeek(moleGran);
        %% Calculate Laplacian molecular thickness
%         thicknessMole = LaplacianThickness(moleGran, 1, 2);
        % volPeek(thicknessMole);colormap jet
        
        %% Save molecular thickness
        niftiwrite(thicknessMole, thicknessMoleNii, ...
            niftiinfo(thicknessFullNii), 'Compressed',true);
        
    end
    
end