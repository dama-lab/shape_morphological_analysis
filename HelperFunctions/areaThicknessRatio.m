function areaThicknessRatio()
    %% calculate surface area to thickness ratio
    
    thicknessNames = {'cortical','granular_WMD','molecular_CSFD'};
    
    %% %%%%%%%%%%%%%%% individual thickness/area
    areaThicknessRatioArray = [];
    
    for l = 1:length(layerTypes)
        layerType = layerTypes{l};
        %% thickness
        thicknessCsv = fullfile(statDir,sprintf('thickness_%s_on_purkinje.csv',thicknessNames{l}));
        thicknessCell = readcell(thicknessCsv);
        thicknessMtx = cell2mat(thicknessCell(3:end,3:end));
        
        %% surface area
        surfaceAreaCsv = fullfile(statDir,sprintf('surface_area_%s_Purkinje.csv',layerType));
        surfaceAreaCell = readcell(surfaceAreaCsv);
        surfaceMtx = cell2mat(surfaceAreaCell(3:end,3:end));
        
        %% thickness/area ratio
        thicknessAreaRatio = thicknessMtx./sqrt(surfaceMtx);
        
        %% save csv 
        csvFname = fullfile(statDir,sprintf('thickness2areaRatio_%s.csv',layerType));
        thicknessAreaRatioCell = thicknessCell;
        thicknessAreaRatioCell(3:end,3:end) = num2cell(thicknessAreaRatio);
        writecell(thicknessAreaRatioCell, csvFname);
        
        %% area/thickness ratio
        areaThicknessRatio = sqrt(surfaceMtx)./thicknessMtx;
        areaThicknessRatioArray = cat(3,areaThicknessRatioArray,areaThicknessRatio);
        
        %% save csv 
        csvFname = fullfile(statDir,sprintf('area2thicknessRatio_%s.csv',layerType));
        areaThicknessRatioCell = thicknessCell;
        areaThicknessRatioCell(3:end,3:end) = num2cell(areaThicknessRatio);
        writecell(areaThicknessRatioCell, csvFname);
        
    end
    
    %% %%%%%%%%%%%%%%% molecular surface area / graular thickness  

    %% Molecular layer surface area
    surfaceAreaCsv  = fullfile(statDir,sprintf('surface_area_%s_Purkinje.csv',layerTypes{3}));
    surfaceAreaCell = readcell(surfaceAreaCsv);
    surfaceMtx = cell2mat(surfaceAreaCell(3:end,3:end));

    %% Granular layer thickness
    thicknessCsv = fullfile(statDir,sprintf('thickness_%s_on_purkinje.csv',thicknessNames{2}));
    thicknessCell = readcell(thicknessCsv);
    thicknessMtx = cell2mat(thicknessCell(3:end,3:end));
    
    %% molecular surface area ratio / graular thickness ratio
    areaThicknessRatio = sqrt(surfaceMtx)./thicknessMtx;
    
    %% save csv 
    csvFname = fullfile(statDir,sprintf('area2thicknessRatio_moleAreaGranThick.csv'));
    areaThicknessRatioCell = thicknessCell;
    areaThicknessRatioCell(3:end,3:end) = num2cell(areaThicknessRatio);
    writecell(areaThicknessRatioCell, csvFname);

    %% granular thickness / molecular surface area ratio ratio
    thicknessAreaRatio = thicknessMtx./sqrt(surfaceMtx);
    
    csvFname = fullfile(statDir,sprintf('thickness2AreaRatio_granThickMoleArea.csv'));
    thicknessAreaRatioCell = thicknessCell;
    thicknessAreaRatioCell(3:end,3:end) = num2cell(thicknessAreaRatio);
    writecell(thicknessAreaRatioCell, csvFname);
    
    %% %%%%%%%%%%%%%%% graular surface area / molecular thickness
    
    %% Granular layer surface area
    surfaceAreaCsv  = fullfile(statDir,sprintf('surface_area_%s_Purkinje.csv',layerTypes{2}));
    surfaceAreaCell = readcell(surfaceAreaCsv);
    surfaceMtx = cell2mat(surfaceAreaCell(3:end,3:end));
        
    %% Molecular layer thickness
    thicknessCsv = fullfile(statDir,sprintf('thickness_%s_on_purkinje.csv',thicknessNames{3}));
    thicknessCell = readcell(thicknessCsv);
    thicknessMtx = cell2mat(thicknessCell(3:end,3:end));
    
    %% molecular surface area ratio / graular thickness ratio
    areaThicknessRatio = sqrt(surfaceMtx)./thicknessMtx;
    
    %% save csv 
    csvFname = fullfile(statDir,sprintf('area2thicknessRatio_granAreaMoleThick.csv'));
    areaThicknessRatioCell = thicknessCell;
    areaThicknessRatioCell(3:end,3:end) = num2cell(areaThicknessRatio);
    writecell(areaThicknessRatioCell, csvFname);
    
    %% graular thickness / molecular surface area ratio ratio
    thicknessAreaRatio = thicknessMtx./sqrt(surfaceMtx);
    
    csvFname = fullfile(statDir,sprintf('thickness2AreaRatio_moleThickGranArea.csv'));
    thicknessAreaRatioCell = thicknessCell;
    thicknessAreaRatioCell(3:end,3:end) = num2cell(thicknessAreaRatio);
    writecell(thicknessAreaRatioCell, csvFname);
    
    
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Zscape for pairwise comparison
    areaThicknessRatios = areaThicknessRatioArray(:,:,2:3);
    %% WT
    areaThicknessRatioWT_gran = areaThicknessRatioArray(1:14,:,2);
    areaThicknessRatioWT_mole = areaThicknessRatioArray(1:14,:,3);

    %% Calculate zscore
    a2tWTboth = [areaThicknessRatioWT_gran,areaThicknessRatioWT_mole]
    a2tMean = mean(a2tWTboth,1);
    a2tstd = std(a2tWTboth,[],1);
    a2tZscore = (a2tWTboth - a2tMean)./a2tstd;
    
    %% figure plot
    figure; imagesc(a2tWTboth)
    
    
    %% TG
    areaThicknessRatioTG_gran = areaThicknessRatioArray(15:end,:,2);
    
        
        