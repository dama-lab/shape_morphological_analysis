function plot_opt = load_opts()
%% image directories
HOME = getenv('HOME');
% Add iso2mesh to the system path
% addpath(genpath(fullfile(HOME,'mial-tools/matlab')));
% addpath(genpath(fullfile(HOME,'Codes/Tools/')));
addpath(genpath(fullfile(HOME,'Codes/faisal-sandbox/students/dma73/projects/cohorts_collection/UCL/mouse_cerebellum')));
addpath(genpath(fullfile(HOME,'Codes/faisal-sandbox/students/dma73/projects/MatlabHelperFunctions')));
rehash % Ref: http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/Installation

%%
rootDir = fullfile(HOME,'Data/Brain_MRI/UCL/Tc1_Cerebellum');
groupMeanDir = fullfile(rootDir,'00_groupwise_average');

maskDir = fullfile(rootDir,'01_mask_groupwise_f3d2');
purkinjeDir = fullfile(rootDir,'07_Purkinje/5_remove_patch');
% statDir = fullfile(rootDir,'08_stats/stats_20191020/remove_missegmented_subjects');
statDir = fullfile(rootDir,'08_stats/stats_20200413');


rawImgDir = fullfile(rootDir,'00_cerebellum_N4_with_mask');
parcellationDir = fullfile(rootDir,'02_parcellation');
sublayerSegDir = fullfile(rootDir,'09_sublayer_parcellation');
sublayerParcellateDir = fullfile(groupMeanDir,'resample_nrr_10/backprop_sublayerParcellation');
granParcDir = fullfile(sublayerSegDir,'/granular_parcellation');
gmRemoveSulciDir = fullfile(rootDir,'03_tissue_label/GM_mask_remove_sulci');
WmGmRemoveSulciDir = fullfile(groupMeanDir,'resample_nrr_10/backprop_Wm1Gm2RemoveSulci');
targetListDir = fullfile(rootDir,'TargetList');

%% Groupwise average
groupMeanNii = fullfile(groupMeanDir,'/mask_manual/average_nonrigid_it_10_nonan.nii.gz');
groupParcellationNii = fullfile(groupMeanDir,'resample_nrr_10/average_nonrigid_it_10_lobuleParcellation.nii.gz');
groupParcellationRecolorNii = fullfile(groupMeanDir,'resample_nrr_10/average_nonrigid_it_10_lobuleParcellationRecolor.nii.gz');
groupPurkinjeNii = fullfile(groupMeanDir,'/segmentations/gmLayerSegKmean_purkinje_lcc_manual.nii.gz');
groupGranularNii = fullfile(groupMeanDir,'/segmentations/gmLayerSegKmean_granularLaplacianManual.nii.gz');
groupGMNii = fullfile(groupMeanDir,'/label_manual/average_nonrigid_it_10_GM_manual.nii.gz');
groupWMNii = fullfile(groupMeanDir,'label_manual/average_nonrigid_it_10_WM_manual.nii.gz');
groupGranularWMNii = fullfile(groupMeanDir,'/segmentations/gmLayerSegKmean_granularPlusWMManual.nii.gz');
groupGranularParcellation = fullfile(groupMeanDir,'/segmentations/gmLayerSegKmean_granular_parcellation.nii.gz');
% groupMeanVol = niftiread(groupMeanNii);
groupCorticalThickness4DNii = fullfile(groupMeanDir,'resample_nrr_10/thicknessCortical/4D.nii.gz');
groupGranularThickness4DNii = fullfile(groupMeanDir,'resample_nrr_10/thicknessGranular/4D.nii.gz');
groupMolecularThickness4DNii = fullfile(groupMeanDir,'resample_nrr_10/thicknessMolecular/4D.nii.gz');
% labelBackPropDir = fullfile(groupMeanDir,'/resample_nrr_10');
granularWmNiiDir = fullfile(groupMeanDir,'resample_nrr_10/backprop_granularWM');

groupRawImgNii = fullfile(groupMeanDir, 'average_nonrigid_it_10_nonan.nii.gz');
groupRawImgCereExtractNii = fullfile(groupMeanDir, 'average_nonrigid_it_10_nonan_cere_extract.nii.gz');
gmLayerSegKmean = fullfile(groupMeanDir,'segmentations/gmLayerSegKmean');
groupGranMoleParcNii = fullfile(groupMeanDir,'segmentations/gmLayerSegKmean_GranMoleParc.nii.gz');

%% Displacement
displaceFildDir = fullfile(groupMeanDir,'resample_nrr_10/displacement');

%% Thickness
thickness_root = fullfile(rootDir,'06_thickness');
thicknessCorticalFullDir = fullfile(thickness_root,'thickness_map');
thicknessCorticalWMDratioDir = fullfile(thickness_root,'thickness_WMD_radio');
thicknessCorticalCSFratioDir = fullfile(thickness_root,'thickness_WMD_radio');
thicknessCortical_CSFD_dir = fullfile(thickness_root,'CSFD');
thicknessCortical_WMD_dir = fullfile(thickness_root,'WMD');
thicknessGranularDir = fullfile(thickness_root,'thickness_granular');
thicknessMoleDir = fullfile(thickness_root,'thickness_molecular');

%% %%%%%%%%%%%%%%%%% Groupwise propagated results
%% Groupwise backpropagated thickness
thicknessDir = fullfile(thickness_root,'groupMeanBackpropagate');
thicknessCorticalFullDir = fullfile(thicknessDir,'thickness_map');
thicknessCortical_CSFD_dir = fullfile(thicknessDir,'CSFD');
thicknessCortical_WMD_dir = fullfile(thicknessDir,'WMD');
%% Groupwise backpropagated sublayer parcellation
corticalParcellateDir = fullfile(groupMeanDir,'resample_nrr_10/backprop_parcellation');
sublayerParcellateDir = fullfile(groupMeanDir,'resample_nrr_10/backprop_sublayerParcellation/');
% sublayerParcellateQCDir = fullfile(groupMeanDir,'resample_nrr_10/backprop_sublayerParcellationQC/');
sublayerParcellateQCDir = fullfile(groupMeanDir,'resample_nrr_10/backprop_sublayerParcellation_inccorrect_colorlabel_but_beautiful/');


%% Surface
surf_dir = fullfile(rootDir,'09_surf');

%% TIV
% TIVCsv = fullfile(rootDir,'/08_stats/stats_20191020/remove_missegmented_subjects/TIV_BV-22.csv');
% TIVCsv = fullfile(rootDir,'/08_stats/stats_20191020/TIV_BV.csv');
TIVCsv = fullfile(rootDir,'/08_stats/stats_20200413/TIV_BV.csv');
[TIV, BV, group, subjectList] = readTIVcsv(TIVCsv);

%% New stat folder for all data
statDirAll = fullfile(rootDir,'08_stats/stats_20200413/');

% result figure dir
resultDir = fullfile(rootDir,'09_figures');
if ~exist(resultDir,'dir'); mkdir(resultDir); end

%% iso2mesh parameters
iso2meshOpt.maxnode = 40000;
iso2meshOpt.radbound = 2;
iso2meshOpt.dofix=1;
iso2meshOpt.isovalues=0.5;
iso2meshOpt.method = 'cgalsurf'; % 'cgalmesh' takes 3x more time (20 sec vs  60 sec)


% create folder for iso2mesh if not already exist
iso2MeshDir = fullfile(surf_dir,'Iso2MeshDir');
if ~exist(iso2MeshDir,'dir'); mkdir(iso2MeshDir); end

%target = '11_27_53-H864-OD-MAS_0';

%% read targetlist
% only 11 WT and 10 TG
% targetListWtFile = fullfile(targetListDir,'TargetList_WT.txt');
% targetListTgFile = fullfile(targetListDir,'TargetList_TG.txt');
% All subjects
targetListWtFile = fullfile(targetListDir,'WT.txt');
targetListTgFile = fullfile(targetListDir,'TG.txt');

targetListWT = readTextLines(targetListWtFile);
targetListTG = readTextLines(targetListTgFile);

% add am empty target to the 10 TG so that it mach with WT number (11)
%targetListTG = cat(1,targetListTG,"");

targetLists = [targetListWT, targetListTG];

groups = {'WT','TG'};
layerTypes = {'full','granular','molecular'};
bgLabel = 0;

%% Surface area / thickness Correlation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Raw thickness
    thicknessFullCsv = fullfile(statDir,'thickness_cortical_on_purkinji.csv');
    thicknessGranCsv = fullfile(statDir,'thickness_granular_WMD_on_purkinje.csv');
    thicknessMoleCsv = fullfile(statDir,'thickness_molecular_CSFD_on_purkinje.csv');
    
    % Raw surface area
    surfaceAreaFullCsv = fullfile(statDir,'surface_area_full_Purkinje.csv');
    surfaceAreaGranCsv = fullfile(statDir,'surface_area_granuler_Purkinje.csv');
    surfaceAreaMoleCsv = fullfile(statDir,'surface_area_molecular_Purkinje.csv');

%% Save all variable in a single .mat file

%%
plot_opt = v2struct;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test new pointCloud class
% https://www.geo.tuwien.ac.at/downloads/pg/pctools/publish/pointCloudIntro/html/pointCloudIntro.html


%% Test extractStruct
% measure + compare: thickness, displacement
% look at measurement on:
%    -Surface displacement (need to scale with TIV?)
%       - whole cortical structure (closed surface)
%       - on Purkinje layer (single layer? dual surface?)
%           need to use non-rigid registration to determine displacement
%   - Thickness
%       - Granular layer
%       - molecular layer

% To-do: extract 3~5 individual significant structure

