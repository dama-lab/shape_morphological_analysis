function fig = cere_surfstat(surfThickMat3D, parcellationVol, surfaceSurfstat, covars, layerTypes, normalizeType, contrast_var, viewAngles)

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
        
        %% Determine normalize type
        if strcmp(normalizeType,'GLM')
            surfThickMat = surfThickMat3D(:,:,layerId);
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

                [slm,stat] = GLMSurfStats(surfThickMat, surfaceSurfstat, covars, contrast_var);
                if strcmp(statMetric,"T")
                    %% [Visualize] T-statistics
                    statsTitle = "T-statistics";
                    [cmap, clim] = statColormap(statMetric);
                    surfaceSurfstat.intensity = slm.t;
                elseif strcmp(statMetric,"P")
                    %% [Visualze] p-value
                    statsTitle = "P-value";
                    % [slm,stat] = GLMSurfStats(surfThickMat./(TIV'), surfaceSurfstat, covars, contrast_var);
                    [cmap, clim, tt, pval] = statColormap(statMetric,stat);
                    surfaceSurfstat.intensity = tt;
                end
                %% actual plot
                if ~exist('viewAngles','var') || isempty(viewAngles)
                    viewAngles = [-87, 54; 95.5, -22.2];
                end
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
    
end