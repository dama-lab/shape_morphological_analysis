
%% 3D volume segmentation
% Ref: https://www.mathworks.com/help/images/3d-volumetric-image-processing.html

addpath(genpath('/home/dma73/Code/projects/MatlabHelperFunctions'));
%%
rootDir = "/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/00_groupwise_average";
rawImgNii = fullfile(rootDir, 'average_nonrigid_it_10_nonan_cere_extract.nii.gz');
wmMaskNii = fullfile(rootDir, 'segmentations/average_nonrigid_it_10_WM_manual.nii.gz');
gmMaskNii = fullfile(rootDir, 'segmentations/average_nonrigid_it_10_GM_manual.nii.gz');
cereMaskNii = fullfile(rootDir, 'mask_manual/average_nonrigid_it_10_nonan.nii.gz');
parcelLabelNii = fullfile(rootDir, 'label_manual/average_nonrigid_it_10_nonan_label_AMBMC-cere-larger-view-GM_sublayer_STEPS_3_8.nii.gz');

gmLayerSegKmean = fullfile(rootDir,'segmentations/gmLayerSegKmean');

rawImg = niftiread(rawImgNii);
wmMask = niftiread(wmMaskNii);
gmMask = niftiread(gmMaskNii);
cereMask = niftiread(cereMaskNii);
parcelLabel = niftiread(parcelLabelNii);
% Remove White Matter
rawImg(wmMask==1) = 0;


%% %%%%%%%%%%%%% Visualilzation volume %%%%%%
%% scrollbar: 
sliceViewer(rawImg);
%% 3-plane https://www.mathworks.com/help/images/ref/orthosliceviewer.html
orthosliceViewer(rawImg);
%% Montage
montage(rawImg,'Indices',1:10:size(rawImg,3)); hold on; 
im = montage(wmMask,'Indices',1:10:size(rawImg,3));
im.AlphaData = 0.5; colormap('jet');
%% Ref: https://www.mathworks.com/help/images/ref/obliqueslice.html
% obliqueslice(rawImg) % not working yet
%% volshow: Ref: https://www.mathworks.com/help/images/ref/volshow.html
% volshow(cropVol(rawImg)); % slow
%% show volume label Ref: https://www.mathworks.com/help/images/ref/labelvolshow.html
powerPC = false; % otherwise will be slow
if powerPC == true
    [rawImgCrop, cropParam] = cropVol(rawImg);
    wmMaskCrop = cropVol(wmMask, cropParam);
    rawImgCropResize = imresize3(rawImgCrop,0.5);
    wmMaskCropResize = imresize3(wmMaskCrop,0.5,'nearest');
    labelvolshow(wmMaskCropResize, rawImgCropResize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K-mean (https://www.mathworks.com/help/images/ref/imsegkmeans3.html)
L = imsegkmeans3(rawImgWhite,5);
L(L==2)=0;
[fig,axesArray] = volPeek(rawImgWhite,L);
figure; orthosliceViewer(L);
%% erode gray matter
SE = strel('sphere',1);
gmMaskErode = imerode(gmMask,SE);
%volPeek(gmMaskErode,uint8(L==4));
%% Purkinje
purkinje = uint8(L==4);
purkinjeGM = uint8(gmMaskErode).*purkinje;
%% Quickcheck
[fig,axesArray] = volPeek(rawImgWhite,purkinjeGM);
figure; orthosliceViewer(purkinjeGM);
%% re-align after resize
relinkAxes(axesArray);
%% Add GM to original K-mean segmentation (Purjinje = 6)
layerSegKmean = purkinjeGM*2 + L;
[fig,axesArray] = volPeek(rawImgWhite,layerSegKmean);
figure; orthosliceViewer(layerSegKmean); colormap jet;
%% Save nifti
niftiwrite(single(layerSegKmean), gmLayerSegKmean, niftiinfo(gmMaskNii), ...
    'Compressed',true);
%% load manual correction
gmLayerSegKmeanManualNii = [char(gmLayerSegKmean),'Manual'];
gmLayerSegKmeanManual = niftiread(gmLayerSegKmeanManualNii);
% sort out labels
gmLayerSegKmeanManual(gmLayerSegKmeanManual==4)=0;
gmLayerSegKmeanManual(gmLayerSegKmeanManual==5)=1;
volPeek(rawImgWhite,gmLayerSegKmeanManual);
niftiwrite(single(gmLayerSegKmeanManual), [gmLayerSegKmeanManualNii,'combined'], niftiinfo(gmMaskNii), ...
    'Compressed',true);
%% Load manual corrected combined EM
gmLayerSegKmeanComnbinedManualNii = [char(gmLayerSegKmean),'CombinedManual.nii.gz'];
gmLayerSegKmeanCombinedManual = niftiread(gmLayerSegKmeanComnbinedManualNii);
volPeek(rawImg,gmLayerSegKmeanCombinedManual);
figure; orthosliceViewer(gmLayerSegKmeanCombinedManual);
%% Purkinje layer
purkinjeSeg = single(gmLayerSegKmeanCombinedManual==6);
volPeek(rawImg,purkinjeSeg);

%% Save manual Purkinje nifti
niftiwrite(single(purkinjeSeg), [char(gmLayerSegKmean),'_purkinje'], niftiinfo(gmMaskNii), ...
    'Compressed',true);

%% reload after manual correction
purkinjeManual = niftiread([char(gmLayerSegKmean),'_purkinje_manual.nii.gz']);
purkinjeManual = single(purkinjeManual==1);
volPeek(rawImg,purkinjeManual);

%% create the GM - Sulci for thickness estimation
purkinjeManual = niftiread('/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/00_groupwise_average/segmentations/gmLayerSegKmean_purkinje_manual.nii.gz');
fissureKmean = uint16(purkinjeManual ==2);
GMremoveSulci = niftiread('/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/00_groupwise_average/resample_nrr_10/average_nonrigid_it_10_GM_mask_remove_sulci.nii.gz');
GMremoveSulci = GMremoveSulci .* (1-fissureKmean);
volPeek(GMremoveSulci); colormap jet

%% save new GMremoveSulci
niftiwrite(single(GMremoveSulci),fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_GM_remove_sulci.nii'),niftiinfo(gmMaskNii),'Compress',true);

%% load manual coorected GM amd plut WM
GMremoveSukciManual = single(niftiread(fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_GM_remove_sulci_manual.nii.gz')));
WM = niftiread(groupWMNii);
WMPlusGMRemoveSulci = WM + 2*GMremoveSukciManual;
volPeek(WMPlusGMRemoveSulci);

%% save WMPlusGMRemoveSulci
niftiwrite(single(WMPlusGMRemoveSulci),fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_WM1_GM2_remove_sulci.nii'),niftiinfo(gmMaskNii),'Compress',true);

%% get Molecular+Sulci
% load files
groupGM = niftiread(groupGMNii);
groupGranular = niftiread(groupGranularNii);
% dilate granular layer by one voxel
SE = strel('sphere',1);
groupGranularDilate = imdilate(groupGranular,SE);
volPeek(rawImg, groupGranularDilate);
%%
groupMolecular = groupGM.*single(1-groupGranularDilate);
groupMolecularExtract = rawImg.*single(groupMolecular);
volPeek(groupMolecularExtract);
%% Kmean-classification-based segmentation
groupMoleKmean = imsegkmeans3(groupMolecularExtract,5);
volPeek(groupMolecularExtract,groupMoleKmean);
%%
figure;orthosliceViewer(groupMoleKmean);colormap jet;

%% Initial sulci
groupSulci = single(groupMoleKmean == 4);
volPeek(permute(groupMolecularExtract,[1,3,2]),permute(groupSulci,[1,3,2]));

%% save sulci
niftiwrite(single(groupSulci),fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_sulci_kmean.nii'),niftiinfo(gmMaskNii),'Compress',true);

%% add sulci to the current manual segmentation
groupGMRemoveSulciManual = niftiread(fullfile(groupMeanDir,'/segmentations/average_nonrigid_it_10_GM_remove_sulci_manual.nii.gz'));
groupSulciCombine = single((single(1-groupGMRemoveSulciManual) + groupSulci) > 0);
volPeek(rawImg,groupSulciCombine);
%% save sulci
niftiwrite(single(groupSulciCombine),fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_sulciPlusBackground.nii'),niftiinfo(gmMaskNii),'Compress',true);
%%groupSulci
volPeek(rawImg,(single(1-groupGMRemoveSulciManual)+groupSulci));

%% save CSF
groupCSF = single(groupMoleKmean == 3);
volPeek(groupMolecularExtract,groupCSF);

%% Add CSF to the manual sulci
groupsulciPlusBackgroundManual = single(niftiread(fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_sulciPlusBackgroundManual.nii.gz')));
groupsulciPlusBackgroundCSF = single((groupCSF + groupsulciPlusBackgroundManual)>0);
volPeek(rawImg,groupsulciPlusBackgroundCSF);
niftiwrite(single(groupsulciPlusBackgroundCSF),fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_sulciPlusBackgroundCSF.nii'),niftiinfo(gmMaskNii),'Compress',true);

%% load manual coorected GM amd plus WM
GMremoveSulciManual = single(niftiread(fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_sulciPlusBackgroundCSFManual.nii.gz')));
WM = niftiread(groupWMNii);
WMPlusGMRemoveSulci = WM + 2*(1-GMremoveSulciManual);
volPeek(WMPlusGMRemoveSulci);

%%
GMremoveSulciManual = single(niftiread(fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_GM_remove_sulci_manual.nii.gz')));
WM = niftiread(groupWMNii);
WMPlusGMRemoveSulci = WM + 2*GMremoveSulciManual;
volPeek(WMPlusGMRemoveSulci);


%% save WMPlusGMRemoveSulci
niftiwrite(single(WMPlusGMRemoveSulci),fullfile(groupMeanDir,'segmentations/average_nonrigid_it_10_WM1_GM2_remove_sulci_manual.nii'),niftiinfo(gmMaskNii),'Compress',true);

%% Largest connected compoinent
purkinjeLCC = LargestConnectedComponent(purkinjeManual);
volPeek(rawImg,purkinjeLCC);

%% save FINAL Purkinje layer segmentation
niftiwrite(single(purkinjeLCC), [char(gmLayerSegKmean),'_purkinje_lcc'], niftiinfo(gmMaskNii), ...
    'Compressed',true);

%% Reload after manual segmentation
purkinjeLCCManual = niftiread([char(gmLayerSegKmean),'_purkinje_lcc_manual.nii.gz']);
volPeek(rawImg, purkinjeLCCManual);
%% %%%%%%%%%%%%%%%%%%%% segment Purkinje/molecular layer
%% Using original layer segmentation (suboptimal)
parcelTwoLayerSeg = parcelLabel;
parcelTwoLayerSeg(parcelTwoLayerSeg==255)=5;
parcelTwoLayerSeg(purkinjeLCCManual==1)=0;
gmThreeLayerSeg = uint16(parcelTwoLayerSeg) + 3 * purkinjeLCCManual;
volPeek(rawImg,gmThreeLayerSeg);
%%
figure;sliceViewer(gmThreeLayerSeg); colormap jet;

%% cleanning thinning Purkinje layer
purkinjeThin = thinning3D(purkinjeLCCManual);
volPeek(purkinjeThin);
% purkinjeLine = bwskel(imcomplement(purkinjeSeg));
% volPeek(uint8(imcomplement(purkinjeSeg)));

%% Prepare cortical labels for Laplacian field
laplacianVol = gmMask - single(purkinjeLCCManual) + 2*wmMask;
volPeek(rawImg,laplacianVol);
%% Calculate Granular layer using Laplacian field
verbose = 1;
[p, granLayer] = LaplaceField(laplacianVol,1,2, verbose);
%%
%granLayer = single(p>0.04 .* p<1);
% granLayer = imfill((p>0.01)*1.* gmMask);
granLayer = imfill(uint8((p>0.01)*1.* gmMask + single(purkinjeLCCManual)>0));
figure; orthosliceViewer(granLayer)
%% save granular layer
niftiwrite(single(granLayer), [char(gmLayerSegKmean),'_granularLaplacian'], niftiinfo(gmMaskNii), ...
    'Compressed',true);
%%
figure;imagesc3D(p)
%% Matlab's Discrete Laplacian operator
% Lap = del2(laplacianVol);
% volPeek(Lap,Lap)

%% Load manually corrected Granular layer
granLayerManual = niftiread([char(gmLayerSegKmean),'_granularLaplacianManual.nii.gz']);
volPeek(rawImg,granLayerManual);

%% add granular layer with white matter
groupWM = niftiread(groupWMNii);
% avoid overlapping voxels
groupWM(logical(uint16(groupWM).*granLayerManual))=0;
groupGranularPlusWM = granLayerManual + uint16(2*groupWM);
volPeek(rawImg,groupGranularPlusWM);

%% save granular + WM
niftiwrite(single(groupGranularPlusWM), [char(gmLayerSegKmean),'_granularPlusWM'], niftiinfo(gmMaskNii), ...
    'Compressed',true);

%% %%%%%%%%%% sublayer parcellation %%%%%%%%%%%%%%%%%
%% Load manual-corrected granular+WM (WM=2, Gran=1)
rawImg = niftiread(groupRawImgNii);
groupGranularPlusWmManual = niftiread([char(gmLayerSegKmean),'_granularPlusWMManual.nii.gz']);
volPeek(rawImg, groupGranularPlusWmManual);

%% Load group GM and group parcellation
groupGM = uint16(niftiread(groupGMNii)); % ad(groupGMNii);
groupParcellation = uint16(niftiread(groupParcellationNii));

%% Create granular + molecular sublayer parcellation
groupGranular = uint16(groupGranularPlusWmManual==1);
groupMolecular = uint16(groupGM) .* uint16(~groupGranularPlusWmManual);
%% incorrect implementation
groupGranMoler = groupGranular + 2*groupMolecular;
groupGranMoleParcellateQC = groupParcellation .* groupGranMoler;
% remove remaining WM dots
groupGranMoleParcellateQC(groupGranMoleParcellateQC==1)=0;
%% save granular + molecular sublayer parcellation parcellation (incorrect but beautiful)
niftiwrite(single(groupGranMoleParcellateQC), [char(gmLayerSegKmean),'_GranMoleParc_QC'], niftiinfo(groupWMNii), ...
    'Compressed',true);

%% Correct implementation
groupGranParcellate = groupGranular .* groupParcellation;
groupMoleParcellate = groupMolecular .* (groupParcellation + 30);
groupGranMoleParcellate = groupGranParcellate + groupMoleParcellate;
% remote elevated background
groupGranMoleParcellate(groupGranMoleParcellate==30)=0;
% remove remaining WM dots
groupGranMoleParcellate(groupGranMoleParcellate==1)=0;
%% groupWmGranMole(groupGM);
volPeek(permute(groupGranMoleParcellate,[1,3,2])); colormap(statColormap("Parcellation",groupGranMoleParcellate));
%%
figure;orthosliceViewer(permute(groupGranMoleParcellate,[1,2,3])); colormap(statColormap("Parcellation",groupGranMoleParcellate));

%% save granular + molecular sublayer parcellation parcellation
niftiwrite(single(groupGranMoleParcellate), [char(gmLayerSegKmean),'_GranMoleParc'], niftiinfo(groupWMNii), ...
    'Compressed',true);

%% Create parcellated gray matter sublayers
groupGM = niftiread(groupGMNii); % ad(groupGMNii);
Group
volPeek(groupGM );

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate gradient
W = gradientweight(purkinjeSeg+gmMask);
volPeek(W);
figure;orthosliceViewer(W);
%% Fast Marching
% seeds = logical(wmMask);
seeds = logical(imdilate(wmMask,SE).*gmMask);
BW = imsegfmm(W,seeds,0.1);
volPeek(purkinjeSeg,uint8(BW));
figure;orthosliceViewer(BW);
%%
gmNoPurkinje = (gmMask - purkinjeSeg)>0;
figure;orthosliceViewer(gmNoPurkinje);

%% Segment granular layer
SE = strel('sphere',1);
SE2 = strel('sphere',2);
figure;
% Initialize granular layer
granularLayer = zeros(size(rawImg));
wmSeed = wmMask;
gmSeed = cereMask==0;
for i = 1:10
    disp(i);
    % dilate for one voxel;
    wmSeed = imdilate(wmSeed,SE);
    gmSeed = imdilate(gmSeed,SE2);
    % remove voxels touching purkinje layer
    wmSeed(logical(wmSeed.*purkinjeSeg))=0;
    gmSeed(logical(gmSeed.*purkinjeSeg))=0;
    % remove voxels touching each other
    wmSeed(logical(wmSeed.*gmSeed))=0;
    gmSeed(logical(wmSeed.*gmSeed))=0;
    % visualize wmSeed evolve progression
    montage(wmSeed,'Indices',50:5:120)
    montage(gmSeed,'Indices',50:5:120)
end
% granularLayer = imdilate(wmMask,SE).*gmNoPurkinje;
%% 
figure;montage(gmSeed);
volPeek(rawImgWhite,wmSeed);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find sulci (https://www.mathworks.com/help/images/ref/imsegfmm.html)
% Compute the weight array based on grayscale intensity differences for mask.
W = graydiffweight(rawImg, logical(gmMask));
volPeek(W);
%% enhance contrast Ref: https://www.mathworks.com/help/matlab/math/powers-and-exponentials.html
W = imadjustn(rawImgWhite.^5);
volPeek(W);
%% volPeek(exp(rawImg));
L = imsegkmeans3(W,5);
L(L==2)=0;
[fig,axesArray] = volPeek(rawImgWhite,L);
orthosliceViewer(L);colormap jet;
%%
purkinje = uint8(L==3);
purkinjeGM = uint8(gmMaskErode).*purkinje;
figure;orthosliceViewer(purkinjeGM);
%%
prctile(rawImg(rawImg~=0),[5,10,25,50,75,90,95])
%%
histogram(rawImg(rawImg~=0))
%%
rawImgWhite = rawImg;
rawImgWhite(rawImgWhite==0) = prctile(rawImg(rawImg~=0),50);
volPeek(rawImgWhite);

%% Purkinje's filter
% Ref: https://www.mathworks.com/matlabcentral/fileexchange/63171-jerman-enhancement-filter
V = planarness3D(rawImgWhite,1:5,[1;1;1],0.75,false);
volPeek(rawImgWhite,V)

%% Frangi's filter
% options.FrangiScaleRange;
V = FrangiFilterPlanar3D(rawImgWhite);
volPeek(rawImgWhite)

%%
% OCTA = '/home/dma73/Data/Retinal_OCT/SankaraNethralaya/RAW_DATA/NII/21/21_angiography.nii.gz';
% II = niftiread(OCTA);
% %%
% options.FrangiScaleRange = [1];
% III = permute(II,[1,3,2]);
% %%
% volPeek(III(:,:,[800:1024]));


%% Fast marching (https://www.mathworks.com/help/images/ref/imsegfmm.html)
%% determine seeds
SE = strel('sphere',1);
seeds = imdilate(wmMask,SE).*gmMask;
% volPeek(seeds);
%% Fast marching
%% speed with raw image
thresh = 0.1;
[BW, D] = imsegfmm(rawImg,logical(seeds),thresh);
volPeek(rawImg,uint8(BW))
%%
% [BW, D] = imsegfmm(W,logical(seeds),thresh);
%% all same speed
thresh = 0.05;
[BW, D] = imsegfmm(ones(size(gmMask)),logical(seeds),thresh);
volPeek(rawImg,uint8(BW))


%%
montage(D>0.05,'Indices',10:5:100)
% volPeek(rawImg,D)
% sliceViewer(D)

%% Superpixel oversegmentation (https://www.mathworks.com/help/images/ref/superpixels3.html)
