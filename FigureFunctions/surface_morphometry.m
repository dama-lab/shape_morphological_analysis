function surface_morphometry()
    %% Plot groupwise purkinje thickness
    % Input:
    %  volMeasure: volume with intensity value of measurements (such as thickness map, parcellation)
    addpath(genpath("/home/dma73/Documents/Code/VHA/vha-dnn/med_deeplearning/cohorts/BrainMRI/Crick/KCL/Dp1Tyb/mouse_cerebellum_UCL_NeuroImage"));


    %% %%%%%%%%%% Preprocessing %%%%%%%%%%%
    
    %% release global variables
    close all;
    plot_opt = load_opts;
    v2struct(plot_opt);
    % load variables
    mask2surfFlag = 1;
    layerTypes = {'full','molecular','granular'};
    thicknesses = {groupCorticalThickness4DNii; ...
                   groupMolecularThickness4DNii; ...
                   groupGranularThickness4DNii};

    %% define function-specific variables
    % matFile to store all inter-media variables (to skip steps and save time)
    matFile = fullfile(surf_dir,'surface_morpormetry_all_28.mat');
    %% create an "empty" mat file if not already exist
    if ~exist(matFile,'file'); save(matFile,'matFile'); end

    %% get the group mean surface layer as the reference targest

    %% [Mask] full GM mask 
    % groupGM = niftiread(groupGMNii);
    % volPeek(groupGranularWM>0);
    
    %% %%%%% generate surface
%     %% [Mask] load granular layer segmentation
%     groupGranular = niftiread(groupGranularNii);
%     %% construct surface with 8 nearest interpolation
%     surface = vol2surface(groupGranular,[],[],[],[],8);
%     figure; surfPlot(surface.vertex,surface.triedge);
    
    %% [Mask] load granular layer + WM segmentation
    groupGranularWM = niftiread(groupGranularWMNii);
    %% construct surface with 8 nearest interpolation
    surface = vol2surface(groupGranularWM,[],[],[],8);
    % figure; surfPlot(surface.vertex,surface.triedge);

    %% load parcellation surface
    % groupParcellation = niftiread(groupParcellationNii);
    groupParcellation = niftiread(groupParcellationRecolorNii);
    parcellationVol = labelRemap(groupParcellation);
    %niftiwrite(parcellationVol,groupParcellationRecolorNii,'Compressed',true);

    %% Vertex parcellation
    vertexParcellation = vol2vertexData(surface.vertex,groupParcellation, [], [], 1);
    % save(matFile,'vertexParcellation','-append');
    %% Setup surface vertex intensity as parcellation
    surface.intensity = vertexParcellation;
    
    %% [QuickCheck] plot parcellation on surface
    % cmap = colormap(colorcube(length(unique(surface.intensity))));
    % cmap = colormap(jet(length(unique(surface.intensity))));
    [cmap,clim,~] = statColormap('Parcellation',parcellationVol);
    surfPlot(surface.vertex, surface.triedge,surface.intensity,[],cmap,clim)
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%
    %% calculate surface thickness for all targets
    surfThickMat2D = [];
    surfThickMat3D = [];
    for layerId = 1:length(layerTypes)
        layerType = layerTypes{layerId};
        fprintf('%s surface thickness:\n', layerType);
        thickness4D = thicknesses{layerId};
        % layer = 'cortical';
        groupThickness4D = niftiread(thickness4D);
        
        %% calculate surface thickness for all targets
        surfThickCell.(layerType) = [];
        tic;
        for targetId = 1:size(groupThickness4D,4)
            target = targetLists{targetId};
            fprintf('    target %s ...\n', target);
            surfThickCell.(layerType)(:,targetId) = vol2vertexData(surface.vertex, groupThickness4D(:,:,:,targetId), [], [], 8);
        end
        toc;
        
        surfThickMat2D = [surfThickMat2D; surfThickCell.(layerType)];
        surfThickMat3D = cat(3,surfThickMat3D, surfThickCell.(layerType));
    end
        
    %% save parameters
    save(matFile,'surface', 'surfThickCell','surfThickMat2D','surfThickMat3D','-append');
    
    %% %%%%%%%%%%% Calculate parcellated mean surface thickness
    %% load thickness 
    load(matFile,'vertexParcellation','surface', 'surfThickCell','surfThickMat2D','surfThickMat3D');
    %% determine labels
    % label remapped:
    %    lob: [1cb  2cb] 3 4/5 6 7 8 9  10 sim Crus1 2  PM Cop PFl F1
    %    raw: [[24 3] 4] 5 6   7 8 9 10 11 12  13    14 15 16  17  18
    %    new: [[20 2] 3] 4 5   6 7 8  9 10 11  12    13 14 15  16  17 (use this)
    % unique(vertexParcellation)
    labels = [20,2:17];

    % Read TIVCsv (to get group name and subject list only)
    [TIVtable, TIVstruct, TIV, BV, group, subjectList] = readTIVcsv(TIVCsv);
    
    % create cell for saving into csv
    structureList = {'1/2Cb','3Cb','4/5Cb','6Cb','7Cb','8Cb','9Cb','10Cb','Sim','Crus 1','Crus 2','PM','Cop','PFl','Fl'};
    
    csvFnames = {'thickness_cortical_on_purkinje.csv'; ...
                 'thickness_molecular_CSFD_on_purkinje.csv'; ...
                 'thickness_granular_WMD_on_purkinje.csv'};
    
    %% Calculate mean thickness for each label, and save to csv
    meanCsv = 0;
    if meanCsv == 1
        for layerId = 1:length(layerTypes)
            %%
            meanThickMtx = [];
            layerType = layerTypes{layerId};
            fprintf('parcellation mean surface thickness of [%s] ...\n', layerType);
            for targetId = 1:length(targetLists(:))
                target = targetLists{targetId};
                % fprintf('    target %s ...\n', target);
                %% combine [1cb (20+2) and 2cb (3)] => 3
                vertexParcellationRemap = vertexParcellation;
                vertexParcellationRemap(ismember(vertexParcellation,[20,2,3])) = 3;
                labelsRemap = 3:17;
                meanThickTarget = parcellatedSurThick(surfThickMat3D(:,targetId,layerId), vertexParcellation, labelsRemap);
                meanThickMtx = [meanThickMtx; meanThickTarget];
            end
            %% combine all cells
            thicknessCell = [['thickness', layerType, structureList]; ...
                             ['groupNo','targets','20+2+3',num2cell(4:17)]; ...
                             [group, num2cell([subjectList,meanThickMtx])]];
            %% Save cell to csv file
            writecell(thicknessCell, fullfile(statDirAll,csvFnames{layerId}));

        end
    end
    %%

    
    %% Compute normalizeed thickness with TIV + group
    surfThickMatNorm3D = [];
    surfThickMatNorm2D = [];
    for layerId = 1:length(layerTypes)
        layerType = layerTypes{layerId};
        fprintf('\n%s:',layerType);

        % Normalize with TIV
        surfThickMat = surfThickMat3D(:,:,layerId);
        % Normalize by (divide) with TIV
        % surfThickMat = surfThickMat./(TIV');
        surfThickMat = GLM(surfThickMat',TIV,group,'wt');
        surfThickMat  = surfThickMat';
        
        %
        surfThickMatNorm2D = [surfThickMatNorm2D; surfThickMat];
        surfThickMatNorm3D = cat(3,surfThickMatNorm3D,surfThickMat);
    end
    %% save normalized surface thickness
    save(matFile,'surfThickMatNorm2D','surfThickMatNorm3D','-append');

    %% Compute normalizeed thickness with TIV only
    surfThickMatNormTIV3D = [];
    surfThickMatNormTIV2D = [];
    for layerId = 1:length(layerTypes)
        layerType = layerTypes{layerId};
        fprintf('\n%s:',layerType);

        % Normalize with TIV
        surfThickMat = surfThickMat3D(:,:,layerId);
        % Normalize by (divide) with TIV
        % surfThickMat = surfThickMat./(TIV');
        surfThickMat = GLM(surfThickMat',TIV,group,'wt');
        surfThickMat  = surfThickMat';
        
        %
        surfThickMatNormTIV2D = [surfThickMatNormTIV2D; surfThickMat];
        surfThickMatNormTIV3D = cat(3,surfThickMatNormTIV3D,surfThickMat);
    end
    %% save normalized surface thickness
    save(matFile,'surfThickMatNormTIV2D','surfThickMatNormTIV3D','-append');
    
    %% %%%%%%%%%% End of preprocessing %%%%%%%%%%%
    
    
    
    
    %% %%%%%%%%%%%%%%%%%%% plot surface thickness difference
    %% statistical analysis / visualization
    v2struct(plot_opt);
    
    clims = [0,20; 0,10;0,10];
    groups = {'WT','TG'};
    % groupIdx = {1:11,12:21}; original 11(WT) + 10(TG) subjects
    groupIdx = {1:14,15:28};
    
    % Load surface
    load(matFile,'surface','surfThickMatNorm3D','surfThickMat3D','surfThickMatNormTIV3D');
    
    normalizeTypes = {'Raw','TIV','GLM'};
    normalizeId = 3;
    normalizeType = normalizeTypes{normalizeId};
    
    %% define covariates
    [~,~,TIV, BV, group, subjectList] = readTIVcsv(TIVCsv);
    covars=[];
    % covars.BV = BV;
    covars.TIV = TIV;
    covars.group = group;
    
    surfaceSurfstat.coord = surface.vertex;
    surfaceSurfstat.tri = surface.triedge;
    % start plot vertex-wise surface morphormetry
    fig = figure; 
    rowNo = 3;
    colNo = 4;
    [ha, pos] = tight_subplot(rowNo+1, colNo, [0,0],[0,0],[0,0]);
    % resize figure
    fig.Position(3)=860;
    fig.Position(4)=800;
    
    rowTitle = {{'Full Cortex'},{'Molecular Layer'},{'Granular Layer'},{'Structural Parcellation'}};
    
    for layerId = 1:length(layerTypes)
        layerType = layerTypes{layerId};
        fprintf('\n%s:',layerType);
        
        %% Normalize with TIV
        if strcmp(normalizeType,'GLM')
            surfThickMat = surfThickMatNorm3D(:,:,layerId);
        elseif strcmp(normalizeType,'TIV')
            surfThickMat = surfThickMatNormTIV3D(:,:,layerId);
        else
            % Normalize by (divide) with TIV
            surfThickMat = surfThickMat3D(:,:,layerId)./(TIV');
        end
        
        %% %%%%%%%% [Statistics] + [Visualize]t-statistics

        %% two camera angles
        %% Statistical analysis
        statsMetrics = {'T','P'};
        for statsId = 1:length(statsMetrics)
            statMetric = statsMetrics{statsId};
            fprintf(' %s-stat..',statMetric);
            
            perspectives = {'front','back'};
            for angleId = 1:length(perspectives)
                perspective = perspectives{angleId};
                fprintf('%s..',perspective);

                if strcmp(statMetric,"T")
                    %% [Visualize] T-statistics
                    statsTitle = "T-statistics";
                    [slm,~] = GLMSurfStats(surfThickMat, surfaceSurfstat, covars);
                    [cmap, clim] = statColormap(statMetric);
                    surfaceSurfstat.intensity = slm.t;
                elseif strcmp(statMetric,"P")
                    %% [Visualze] p-value
                    statsTitle = "P-value";
                    [slm,stat] = GLMSurfStats(surfThickMat./(TIV'), surfaceSurfstat, covars);
                    [cmap, clim, tt, pval] = statColormap(statMetric,stat);
                    surfaceSurfstat.intensity = tt;
                end
                %% actual plot
                viewAngles = [-87, 54; 95.5, -22.2];
                axisId = (layerId-1)*colNo + (statsId-1)*length(perspectives)+angleId;
                axes(ha(axisId));
                if strcmp(perspective,"front")
                    viewAngle = viewAngles(1,:); % [-87, 54];
                elseif strcmp(perspective,"back")
                    viewAngle = viewAngles(2,:); % [95.5, -22.2];
                end
                
                [surfFig,TR, colormaps] = surfPlot(...
                    surfaceSurfstat.coord, surfaceSurfstat.tri, surfaceSurfstat.intensity, [], ...
                    cmap, clim, viewAngle,[]);
                
                %% add figure marks
                %% title for statsitical analysis
                % title.Position(1): height
                % title.Position(2): distance towards left
                if ismember(axisId,[2,4])
                    t = title(gca, sprintf('%s',statsTitle));
                    t.Position(1) = 50;
                    t.Position(2) = 0;
                    t.FontAngle = 'italic';
                    t.FontSize = 10;
                end
                %% title for layer type
                if ismember(axisId,[3,7,11])
                    t = title(gca, rowTitle{layerId});
                    t.Position(1) = 130;
                    t.Position(2) = 265;
                    t.FontSize = 11;
                end
                % title(layerType);
                %colorbar;
                % caxis([-4,4]); colormap(cmap);
                % Set nan values to transparent: 
                
                % [Visualze] q-value (False discovery rate)
                % figure; SurfStatView(stat.qval, surfaceSurfstat, 'WT v.s. TG removing ICV');

                % Inflation
                % SurfStatInflate(surfaceSurfstat)
            end
        end
    end
    
    %% update clim for the T-map
    for haId = 1:length(ha)
        if ismember(haId,[1,2,5,6,9,10])
            climT = [-6,6];
            caxis(ha(haId),climT);
        end
    end
    
    %% colorbar for T
    % 1: horizontal location (related to the whole figure)
    % 2: vertical location  (related to the whole figure)
    % 3: width
    % 4: height
    % [colormaps, clim] = statColormap('T');
    % axes(ha(10));
    % colormap(gca,colormaps);
    cbt = colorbar(ha(10),'south','Ticks',[min(climT),min(climT)/6,max(climT)/6,max(climT)],'Box','on');
    cbt.Position(1)=0.13;
    cbt.Position(2)=0.28;
    cbt.Position(4)=0.018;
    cbt.FontSize = 8;
    cbt.TickLength = 0;
    set(cbt,'xAxisLocation','bottom');
    
    %% colorbar for P
    axes(ha(12));
    axisP = gca;
    clim = axisP.CLim;
    cbp = colorbar(ha(12),'south');
    cbp.Position(1) = 0.6;
    cbp.Position(2)=0.28;
    cbp.Position(4)=0.018;
    cbp.Ticks = [clim(1),6,7.5,clim(2)];
    cbp.TickLabels = [pval.thresh, 0, pval.thresh, 0];
    cbp.TickLength = 0;
    set(cbp,'xAxisLocation','bottom');
    
    %% plot parcellation
    surfaceSurfstat.intensity = vol2vertexData(surfaceSurfstat.coord,parcellationVol, [], [], 1);
    [cmap,clim,~] = statColormap('Parcellation',parcellationVol);
    axes(ha(13));
    surfPlot(surfaceSurfstat.coord, surfaceSurfstat.tri, surfaceSurfstat.intensity, [], ...
                    cmap, clim, viewAngles(1,:));
    % turn shading interp off
    axes(ha(14));
    opts.shading='flat';
    surfPlot(surfaceSurfstat.coord, surfaceSurfstat.tri, surfaceSurfstat.intensity, [], ...
                    cmap, clim, viewAngles(2,:), opts);
                
    

    %% Add parcellation Legend
    axes(ha(15));
    labelNo=16;
    cmap16 = cmap(3:2+labelNo,:);
    x = ones(1,labelNo); %[zeros(1,labelNo);ones(1,labelNo)];
    y = x; %2*(1:labelNo);
    for p = 1:size(x,2)
        plot(x(p),y(p),'.','MarkerSize',30,'Color',cmap16(p,:));
        hold on;
    end
    labelName={'1Cb','2Cb','3Cb','4/5Cb','6Cb','7Cb','8Cb','9Cb','10Cb',...
        'Sim','Crus 1','Crus 2','PM','Cop','PFl','Fl'};
    lgd = legend(labelName,'NumColumns',4,'FontSize',10,'Location','west');
    lgd.Position(1) = 0.55 ;
    ha(15).Visible = 'off';   
    
    %% remove the last axis
    ha(16).Visible = 'off';
    
    %% Save figure
    figFname = fullfile(resultDir,'/surface_stat_all_28.png');
    saveMethod = "export_fig"; % "native"
    if strcmp(saveMethod, "export_fig")
        % preferred way
        export_fig(figFname, '-m2');
    elseif strcmp(saveMethod, "native")
        % Native way to save figure
        MatlabVersion = version;
        MatlabVersion =MatlabVersion(end-5:end-2);
        if str2double(MatlabVersion) >= 2020
            t.Units = 'inches';
            t.OuterPosition = [0 0 10 25];
            exportgraphics (fig, [figFname,'.png'],'Resolution', 200);
        else
            print(fig, figFname, '-dpng', '-r200');
        end
    end
    
    %% stop here
    return;
    
    %% %%%%%%%%%%%%%%%%% Plot mean surface 
    surfaceThicknesses=[];
    fig = figure; 
%     [ha, pos] = tight_subplot(2,3,[0,0],[0,0],[0,0]);
     % tiledlayout(3,2);
    for layerId = 1:length(thicknesses)
        thickness4D = thicknesses{layerId};
        % layer = 'cortical';
        groupThickness4D = niftiread(thickness4D);
        %%
        for groupId = 1:length(groups)
            %% extract each group
            groupThickness = groupThickness4D(:,:,:,groupIdx{groupId});
            %% calculate the mean
            groupThicknessMean = mean(groupThickness,4);
            %% [To-Do] (Normalize with TIV)
            %% construc/plot surface for each group
            surface.thickness = vol2vertexData(surface.vertex, groupThicknessMean, [], [], 8);
            surfaceThicknesses = [surfaceThicknesses, surface.thickness];
            subplot(2,3,(groupId-1)*length(thicknesses)+layerId);
            %ha((length(thicknesses)-1)*groupId+thicknessId); % nexttile;
            clim = clims(layerId,:);
            surfPlot(surface.vertex, surface.triedge, surface.thickness,[],[],clim);
        end
    end

    
    %% %%%%%%%%%%%%%%%%%%%%%%% Distance <= Displacement %%%%%%%%%%%%%%%%%%%
    %% Calculate the surface distance from displacement field for all targets
    surfDistMat2D = [];
    surfDistMat3D = [];
    % layerTypes = {'full','granular','molecular'};
    %%
    for dispId = 2:length(layerTypes)
        layerType = layerTypes{dispId};
        fprintf('%s surface displacement:\n', layerType);
        %% 
        surfDistCell.(layerType) = [];
        for targetId = 1:(length(targetLists)-1)
            %%
            target = targetLists(targetId);
            fprintf('    target %s ...\n', target);
            volDisplacement = squeeze(niftiread(fullfile(displaceFildDir,target)));
            % calculate surface displacement => distance
            tic;
            surfDistance = vol2surfDisplacement(surface.vertex, surface.triedge, volDisplacement);
            toc;
            surfDistCell.(layerType)(:,targetId) = surfDistance;
            surfPlot(surface.vertex,surface.triedge,surfDistance,[],[],[-1,1])
        end
        %%
        surfDistMat2D = [surfDistMat2D, surfDistCell.(layerType)];
        surfDistMat3D = cat(3, surfDistMat3D, surfDistCell.(layerType));
    end
    
%         %% extract each group
%         groupThicknessWT = groupThickness4D(:,:,:,1:11);
%         groupThicknessTG = groupThickness4D(:,:,:,12:21);
%         %% calculate the mean
%         groupThicknessWTMean = mean(groupThicknessWT,4);
%         groupThicknessTGMean = mean(groupThicknessTG,4);
%         %%
%         % niftiwrite(groupThicknessWTMean,fullfile(groupMeanDir,'/resample_nrr_10/corticalThickness/meanWT','Compressed',true));
%         %%
%         % volPeek(groupThicknessWTMean);colormap jet
%         %% [To-Do] (Normalize with TIV)
%         surfThickWT = vol2vertexData(surface.vertex,groupThicknessWTMean, [], [], 8);
%         surfThickTG = vol2vertexData(surface.vertex,groupThicknessTGMean, [], [], 8);
%         %% construc/plot surface for each group 
%         %% WT
%         surfThicknessWT = vol2vertexData(surface.vertex,groupThicknessWTMean, [], [], 8);
% %         intStat = statsMetrics(surface.intensity);
% %         clim = [max(intStat.min,intStat.mean-3*intStat.std), ...
% %                 min(intStat.max,intStat.mean+3*intStat.std)];
%         subplot(1,2,1); surfPlot(surface.vertex, surface.triedge, intensityWT,[],[],clim);
%         %% TG
%         surfThicknessTG = vol2vertexData(surface.vertex,groupThicknessTGMean, [], [], 8);
%         subplot(1,2,2); surfPlot(surface.vertex, surface.triedge, surfThicknessTG, [],[],clim);
%     
%     
%     %%
%     figure;pcshow(pointCloud(surface.vertex,'Intensity',surface.intensity), ...
%         'MarkerSize',1000);
    
%     %% %%%%%%%%% load structural parcellation
%     groupParcellation = uint16(niftiread(groupParcellationNii));
%     
%     %% parcellate the granular
%     groupGranularParcellation = uint16(groupGranularWM>0) .* groupParcellation;
%     volPeek(permute(groupGranularParcellation,[1,3,2]));colormap jet
%     
%     %% extract specific lobule
%     lobuleNo = 9;
%     groupGranularLobule = (groupGranularParcellation == lobuleNo);
%     volPeek(permute(groupGranularLobule,[1,3,2])); colormap jet
%     groupGranularSurf = vol2surface(groupGranularParcellation==lobuleNo, 0, '', 0);
    
    %% convert vol2surface
    
    %%
    figure;
    % pcshow(groupGranularSurf.vertex)       
    % surfPlot(groupGranularSurf.vertex,groupGranularSurf.triedge)
     [surfplot, TR] = surfPlotVol(groupGranularWM>0);
    
    
    %% %%%%%%%%%%%%%%%%%%% convert surface into mesh
    %% convert to PointCloud
    pc = pointCloud(groupGranularSurf.vertex);
    %% pcdownsample
    pcDownsample = pcdownsample(pc, ...
        'gridAverage', 2);
    groupGranularSurf.downsampleVertex = pcDownsample.Location;
    pcshow(groupGranularSurf.downsampleVertex)
    %% Method one (need to be smoothed): MyCrust (Ref: https://www.mathworks.com/matlabcentral/profile/authors/1548775-luigi-giaccari)
    % close version:
    % groupPurkinjeSurf.triedge = MyRobustCrust(groupPurkinjeSurf.vertex);
    % groupPurkinjeSurf.downsampleTriedge = MyRobustCrust(groupPurkinjeSurf.downsampleVertex);
    
    % Open version:
    groupGranularSurf.triedge = MyCrustOpen(groupGranularSurf.vertex);
    % groupPurkinjeSurf.downsampleTriedge = MyCrustOpen(groupPurkinjeSurf.downsampleVertex);
    
    %% fix the mesh
    % vertexfix = groupPurkinjeSurf.downsampleVertex;
    % triedgefix = groupPurkinjeSurf.downsampleTriedge;
    vertexfix = groupGranularSurf.vertex;
    triedgefix = groupGranularSurf.triedge;
    
    [vertexfix,triedgefix] = meshcheckrepair(vertexfix,triedgefix,'dup');
    [vertexfix,triedgefix] = meshcheckrepair(vertexfix,triedgefix,'deep');
    [vertexfix,triedgefix] = meshcheckrepair(vertexfix,triedgefix,'isolated');
    %
    surfPlot(vertexfix,triedgefix)
    
    %% find the vertexnorm of the initial surface
    vnormal = vertexNormal(triangulation(double(triedgefix),vertexfix));
    
    %% Method 2 IBSL3_3DT/IBSL3_3DTRI (Ref: https://www.mathworks.com/matlabcentral/profile/authors/3793616-mohammad-rouhani)
    % Normalize
    vertexNorm = normalizeCoordinate(vertexfix);
    % trancated:
    L = 10; %regularization parameter; increase it for a coarser surface.
    Size = 100; %IBS size can be increased to 100.
    P = IBSL3_3D_T(0.01,Size,L,vertexNorm, 0, vnormal);
    % P = IBSL3_3D_T(0.01,Size,L,S,0,N);
    % figure;surfPlot(S,P+1);
    %%
    figure;surfPlot(S,T+1)
    % figurelIBSLevelSurf
    
    %% 
    figure; surfPlot(groupGranularSurf.vertex,t)
    
    %% save group surface
    save(matFile,'groupPurkinjeSurf','-append');

    return;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Extract individual purkinje layer and convert to point cloud
    data_dir = purkinjeDir;
    dirlist = dir(data_dir);
    targetlist = dirlist(3:end);


    %%
    for target_no = 1:length(targetlist)
        %% Register point cloud to the reference point cloud (get the transformation)
        target_no = 1;
        target_id = targetlist(target_no);
        %% Resample using the transformation

        %% 
    end
end