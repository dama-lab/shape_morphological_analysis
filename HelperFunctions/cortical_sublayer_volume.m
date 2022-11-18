function cortical_sublayer_volume(plot_opt)
    %% load all variables
    v2struct(plot_opt);
    %% Target List
    targetNo = length(targetLists(:));
    
    %% Read TIVCsv
    [TIV, BV, group, subjectList] = readTIVcsv(TIVCsv);
    
    %% create cell for saving into csv
    structureList = {'1/2Cb','3Cb','4/5Cb','6Cb','7Cb','8Cb','9Cb','10Cb','Sim','Crus 1','Crus 2','PM','Cop','PFl','Fl'};
    
    %% %%%%%%%%%%%%%%%%%%%%
    %% whole cerebellar volume
    volCereArray = [];
    for targetId = 1:targetNo
        target = targetLists{targetId};
        fprintf('%s\n',target);
        %% load parcellation labels
        segVolNii = fullfile(corticalParcellateDir,target);
        segVol = niftiread(segVolNii);
        niftiheader = niftiinfo(segVolNii);
        voxeldim = niftiheader.PixelDimensions;
        volTarget = extractVol(segVol,[],voxeldim);
        % Calculate cerebellar volume (exclude background voxels)
        volCere = sum(volTarget(2:end));
        volCereArray = [volCereArray; volCere];
    end
    %%
    writematrix(volCereArray, fullfile(statDirAll,'volume_cerebellar_mask.csv'));
    
        
    
    %% %%%%%%%%%%%%%%%%%%%%
    %% Full cortical parcellation
    volMtx = [];
    for targetId = 1:targetNo
        target = targetLists{targetId};
        fprintf('%s\n',target);
        %% load parcellation labels
        segVolNii = fullfile(corticalParcellateDir,target);
        segVol = niftiread(segVolNii);
        %% label remapped:
        % lob: [1cb  2cb] 3 4/5 6 7 8 9  10 sim Crus1 2  PM Cop PFl F1
        % raw: [[24 3] 4] 5 6   7 8 9 10 11 12  13    14 15 16  17  18
        % new: [[20 2] 3] 4 5   6 7 8  9 10 11  12    13 14 15  16  17 (use this)
        volTarget = extractVol(segVol, [20,2:17]);
        volMtx = [volMtx; volTarget];
    end
    
    %% %%%%%%%%%%%% save full parcellation
    
    % combine [1cb (20+2) and 2cb (3)]
    volMtcRemap = [sum(volMtx(:,1:3),2), volMtx(:,4:end)];
        
    % combine all cells
    volCell = [['volume','full_cortical',structureList]; ...
               ['voxelNo','targets','20+2+3',num2cell(4:17)]; ...
               [group, num2cell([subjectList,volMtcRemap])]];
    %% write cell to csv
    writecell(volCell, fullfile(statDirAll,'volume_full_raw.csv'));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Sublayer cortical parcellation
    granMtx = [];
    moleMtx = [];
    
    for targetId = 1:targetNo
        target = targetLists{targetId};
        fprintf('\n%s ',target);
        %% load parcellation labels
        segVolNii = fullfile(sublayerParcellateDir,target);
        segVol = niftiread(segVolNii);
        %% label remapped:
        % lob: [1cb  2cb] 3 4/5 6 7 8 9  10 sim Crus1 2  PM Cop PFl F1
        % raw: [[24 3] 4] 5 6   7 8 9 10 11 12  13    14 15 16  17  18 (granular)
        %     [[54,33]34]35 ... (Molecular layer)
        % new: [[20 2] 3] 4 5   6 7 8  9 10 11  12    13 14 15  16  17
        %% granular layer
        fprintf('granular ');
        granTarget = extractVol(segVol, [24,2:17]);
        granMtx = [granMtx; granTarget];
        %% molecular layer
        fprintf('molecular ');
        moleTarget = extractVol(segVol, 30+[24,2:17]);
        moleMtx = [moleMtx; moleTarget];
    end
    
    %% %%%%%%%%%%% save Granular layer
    % combine [1cb (20+2) and 2cb (3)]
    granMtcRemap = [sum(granMtx(:,1:3),2), granMtx(:,4:end)];
        
    % combine all cells
    granCell = [['volume','granular_layer',structureList]; ...
               ['voxelNo','targets','24+3+4',num2cell(5:18)]; ...
               [group, num2cell([subjectList,granMtcRemap])]];

    %% write cell to csv
    writecell(granCell, fullfile(statDirAll,'volume_granular_parcellation.csv'));

    %% %%%%%%%%%%% save Molecular layer
    % combine [1cb (20+2) and 2cb (3)]
    moleMtcRemap = [sum(moleMtx(:,1:3),2), moleMtx(:,4:end)];
        
    % combine all cells
    moleCell = [['volume','granular_layer',structureList]; ...
               ['voxelNo','targets','24+3+4',num2cell(5:18)]; ...
               [group, num2cell([subjectList,moleMtcRemap])]];

    %% write cell to csv
    writecell(moleCell, fullfile(statDirAll,'volume_molecular_parcellation.csv'));

    
    
    
    
    
    
    
    
    
    