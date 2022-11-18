function groupwise_purkinje_thickness(plot_opt)
%% Plot groupwise purkinje thickness

% Convert Purkinje mask to surface plot
mask2surfFlag = 1;

labelNo = 9;

%% release global variables
v2struct(plot_opt);

%% define function-specific variables
% matFile to store all inter-media variables (to skip steps and save time)
matFile = fullfile(surf_dir,'purkinje_thickness.mat');
% greate an "empty" mat file if not already exist
if ~exist(matFile,'var'); save(matFile,'matFile'); end

%% get the group mean purkinje layer as the reference targest
if mask2surfFlag == 1    
    %% Full purkinje layer
    groupPurkinje = niftiread(groupPurkinjeNii);
    % volPeek(groupPurkinje==1);
    
    
    %% thinining the Purkinje layer
    groupPurkinjeThin = thinning3D(groupPurkinje);
    % volPeek(permute(groupPurkinjeThin,[1,3,2]));
    
    %% load structural parcellation
    groupParcellation = niftiread(groupParcellationNii);
    
    %% parcellate the Purkinje
    groupPurkinjeParcellation = groupPurkinjeThin .* groupParcellation;
    % volPeek(permute(groupPurkinjeParcellation,[1,3,2]));colormap jet
    
    %% extract specific lobule
    lobuleNo = 9;
    groupPurkinjeLobule = (groupPurkinjeParcellation == lobuleNo);
    volPeek(groupPurkinjeLobule); colormap jet
    
    %% convert vol2surface
    groupPurkinjeSurf = vol2surface(groupPurkinjeParcellation==lobuleNo, 0, '', 1);
    figure;pcshow(groupPurkinjeSurf.vertex)       
        
    
    %% %%%%%%%%%%%%%%%%%%% convert surface into mesh
    %% convert to PointCloud
    pc = pointCloud(groupPurkinjeSurf.vertex);
    %% pcdownsample
    pcDownsample = pcdownsample(pc, ...
        'gridAverage', 2);
    groupPurkinjeSurf.downsampleVertex = pcDownsample.Location;
    pcshow(groupPurkinjeSurf.downsampleVertex)
    %% Method one (need to be smoothed): MyCrust (Ref: https://www.mathworks.com/matlabcentral/profile/authors/1548775-luigi-giaccari)
    % close version:
    % groupPurkinjeSurf.triedge = MyRobustCrust(groupPurkinjeSurf.vertex);
    % groupPurkinjeSurf.downsampleTriedge = MyRobustCrust(groupPurkinjeSurf.downsampleVertex);
    
    % Open version:
    groupPurkinjeSurf.triedge = MyCrustOpen(groupPurkinjeSurf.vertex);
    % groupPurkinjeSurf.downsampleTriedge = MyCrustOpen(groupPurkinjeSurf.downsampleVertex);
    
    %% fix the mesh
    % vertexfix = groupPurkinjeSurf.downsampleVertex;
    % triedgefix = groupPurkinjeSurf.downsampleTriedge;
    vertexfix = groupPurkinjeSurf.vertex;
    triedgefix = groupPurkinjeSurf.triedge;
    
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
    figure; surfPlot(groupPurkinjeSurf.vertex,t)
    
    %% save group surface
    save(matFile,'groupPurkinjeSurf','-append');
end

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
