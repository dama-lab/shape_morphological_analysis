
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
noRow = 6;
noCol = 8;

% get subpIdx
subpIdx = transpose(reshape(1:noRow*noCol, [noCol,noRow]));

%% get [groupTargetIdx] matrix
%% WT
wtIdx = [1:6;1:6;7:12;7:12]';
wtIdx(wtIdx==12)=0;
% TG
tgIdx = [[1:5,0];[1:5,0];[6:10,0];[6:10,0]]';
% combined
targetIdx = [wtIdx,tgIdx];

% groupIdx combined
groupIdx = [ones(6,4), 2*ones(6,4)];
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
    target = targetLists(targetId,groupId);
    
    %%
    if mod(sbpltId,2)
        fprintf('%s(%d) - %s\n', groups{groupId},targetId,target);
        
        %% load parcellation images
        % parcellateNii = fullfile(parcellationDir,target);
        parcellateNii = fullfile(sublayerParcellateDir,target);
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
        padding = 5;
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
    %rawImgMidSlice = squeeze(rawCrop(:, :, round(size(rawCrop,3)/2+sliceShift)));
    imagesc(ax, rawImgMidSlice);
    colormap(ax, gray);
    axis(ax, 'off'); % , 'vis3d'
    
    %% plot overlay plotting if subplot number = odd
    if ~mod(sbpltId,2)
        parcellateMidSlice = squeeze(parcellateCrop(:, round(size(parcellateCrop,2)/2)+sliceShift, :));
        %parcellateMidSlice = squeeze(parcellateCrop(:,:,round(size(parcellateCrop,3)/2+sliceShift)));
        
        %plotOverlay(parcellateMidSlice, ax);
        ax1 = axes;
        plt = imagesc(ax1,parcellateMidSlice);
        axis(ax1,'off')
        plt.AlphaData = ~isnan(parcellateMidSlice);
        alpha=0.1;
        plt.AlphaData(plt.AlphaData==1)=alpha;
        colormap(ax1,colormapFull);
        linkprop([ax,ax1],'Position');
        % title(ax1,target);
    end
    
    switch sbpltId
        case 3
            tt = title(ha(sbpltId),'WildType');
            tt.Position(1) = 0;
            tt.FontSize = 12;
        case 7
            tt = title(ha(sbpltId),'Transcromosomic');
            tt.Position(1)=0;
            tt.FontSize = 12;
    end
    
end

% resize figure (keep hight/width ratio)
whRatio = fig.Position(4)/fig.Position(3);
fig.Position(3) = 1000;
fig.Position(4) = fig.Position(3) * whRatio;

%% save figure 
figFname = fullfile(resultDir,'segmentation_qc.png');
save_figure(fig,figFname);
