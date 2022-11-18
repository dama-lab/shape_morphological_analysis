function individualSegmentationQC_all_28_subjects()
%% Load structures
plot_opt = load_opts();
v2struct(plot_opt);

%% Generate colormaps
parcellationVol = niftiread(fullfile(parcellationDir,targetLists{1}));
colormaps = statColormap("Parcellation",parcellationVol);
colormapReverse = [colormaps(:,3) colormaps(:,2) colormaps(:,1)];
colormapPermute = [colormaps(:,2) colormaps(:,3) colormaps(:,1)];
% colormapShuffle = [colormaps(:,1) colormaps(:,3) colormaps(:,2)];
colormapReverse(10,:) = colormapPermute(2,:);
colormapFull = [colormaps; colormapReverse(1:9,:)];

%% Mapping the groupwise target index to the subplot index
noRow = 7;
noCol = 8;

% get subpIdx
subpIdx = transpose(reshape(1:noRow*noCol, [noCol,noRow]));

%% get [groupTargetIdx] matrix
% WT
wtIdx = [1:7;1:7;8:14;8:14]';
% TG
tgIdx = wtIdx;
% combined
targetIdx = [wtIdx,tgIdx];

% groupIdx combined
groupIdx = [ones(7,4), 2*ones(7,4)];
groupTargetIdx = cat(3, groupIdx, targetIdx);

%% Initialize figure
fig = figure;
[ha,] = tight_subplot(noRow,noCol,[0,0],[0,0.05],[0,0]);
%
groups = {'wt','tg'};
for sbpltId = 1:length(ha)
    %% determing group targetId
    groupId = groupIdx(subpIdx==sbpltId);
    targetId = targetIdx(subpIdx==sbpltId);
    %% SKipping non-existing targets
    if targetId == 0
        ha(sbpltId).Visible='off';
        continue; 
    end
    %% determine targets
    target = targetLists{targetId,groupId};
    
    %%
    if mod(sbpltId,2)
        fprintf('%s(%d) - %s\n', groups{groupId},targetId,target);
        
        %% load parcellation images
        % parcellateNii = fullfile(parcellationDir,target);
        % parcellateNii = fullfile(sublayerParcellateQCDir,target);
        % parcellateNii = fullfile(sublayerParcellateDir,target);
        parcellateNii = fullfile('/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/00_groupwise_average/resample_nrr_10/backprop_sublayerParcellationQC/',target);
        % convert to double so that can assign 0 to nan
        parcellateVol = double(niftiread(parcellateNii));
%         % remove WH
%         parcellateVol(parcellateVol==1)=0;
%         % remove boundary labels
%         parcellateVol(parcellateVol==2)=22;
        % convert backgroun 0 as nan
        bgLabel = 0;
        parcellateVol(parcellateVol == bgLabel) = NaN;

        %% crop segmentation
        padding = 40;
        [parcellateCrop, cropParam] = cropVol(parcellateVol, [], padding);

        %% load raw vol
        rawNii = fullfile(rawImgDir,target);
        rawVol = niftiread(rawNii);

        %% crop/pad raw vol
        rawCrop = cropVol(rawVol, cropParam);
    end
    
    sliceShift = 0;
    %% plot raw image
    ax = ha(sbpltId);
    rawImgMidSlice = squeeze(rawCrop(:, round(size(rawCrop,2)/2)+sliceShift, :));
    rawImgMidSlice = rot90(permute(rawImgMidSlice,[2,1,3]),2);
    %rawImgMidSlice = squeeze(rawCrop(:, :, round(size(rawCrop,3)/2+sliceShift)));
    imagesc(ax, rawImgMidSlice);
    colormap(ax, gray);
    axis(ax, 'vis3d', 'off'); % , 'vis3d'
    % daspect(ax,'auto');
    % ax.DataAspectRatio = 
    
    %% plot overlay plotting if subplot number = odd
    if ~mod(sbpltId,2)
        parcellateMidSlice = squeeze(parcellateCrop(:, round(size(parcellateCrop,2)/2)+sliceShift, :));
        %parcellateMidSlice = squeeze(parcellateCrop(:,:,round(size(parcellateCrop,3)/2+sliceShift)));
        parcellateMidSlice = rot90(permute(parcellateMidSlice,[2,1,3]),2);
        %plotOverlay(parcellateMidSlice, ax);
        ax1 = axes;
        plt = imagesc(ax1,parcellateMidSlice);
        axis(ax1,'vis3d','off');
        plt.AlphaData = ~isnan(parcellateMidSlice);
        alpha=0.1;
        plt.AlphaData(plt.AlphaData==1)=alpha;
        colormap(ax1,colormapFull);
        caxis(ax1,[0,58]);
        linkprop([ax,ax1],'Position');
        % title(ax1,target);
    end
    
    switch sbpltId
        case 3
            tt = title(ha(sbpltId),'WildType');
            tt.Position(1) = 0;
            tt.FontSize = 12;
        case 7
            tt = title(ha(sbpltId),'Tc1');
            tt.Position(1)=0;
            tt.FontSize = 12;
    end
    
end

%% resize figure (keep hight/width ratio)
whRatio = fig.Position(4)/fig.Position(3);
fig.Position(3) = 840;
fig.Position(4) = 760; %fig.Position(3) * whRatio;

%% save figure 
figFname = fullfile(resultDir,'segmentation_qc_all_28_original_aspect_ratio.png');
save_figure(fig,figFname);
