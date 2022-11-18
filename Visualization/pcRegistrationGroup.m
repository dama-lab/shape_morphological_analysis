function tformGroup = pcRegistrationGroup(surfGroup,refTarget)
%% Groupwise Point-Cloud Registration
%    surfGroup: containers.Map
%    [Option 1]refSurf: reference surface (Recommented: from groupwise registration average)
%    [Option 2]refTarget: target as reference surface

targetList = surfGroup.keys;

%% check if refTarget exist, conver to point cloud
if ~isKey(surfGroup,refTarget)
    error('did not find reference surf');
else
    fprintf('using %s as reference image\n', refTarget)
    surfRef = surfGroup(refTarget);
    % convert to point cloud
    ptCloudRef = pointCloud(surfRef);
end

%% Registration to the reference images
tformGroup = containers.Map;
distGroup = containers.Map;
for targetId = 1:length(targetList)
    target = targetList{targetId};
    fprintf('Registering: %s\n',target);
    %% skip reference target
    if strcmp(target,refTarget); continue; end
    
    %% convert surf to point cloud
    surfCurrent = surfGroup(target);
    ptCloudCurrent = pointCloud(surfCurrent);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % point cloud registration. Ref: https://www.mathworks.com/help/vision/examples/3-d-point-cloud-registration-and-stitching.html
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% downsample to speed up and improve the accuracy of the registration
    gridSize = 0.1;
    fixed  = pcdownsample(ptCloudRef,'gridAverage', gridSize);
    moving = pcdownsample(ptCloudCurrent,'gridAverage', gridSize);
    
    %% rigid registration (ICP)
    tformRigid = pcregistericp(moving, fixed, 'Metric', 'pointToPoint', 'Extrapolate', true);
    % ptCloudRigid = pctransform(moving, tformRigid);
    
    tformGroup(target) = tformRigid;
    
    
%     %% affine registration (CPD)
%     tformAffine = pcregistercpd(ptCloudRigid,fixed,'Transform','Affine','Verbose',true);
%     ptCloudAffine = pctransform(ptCloudRigid, tformAffine);
%     
%     %% affine registration (CPD, no rigid initialization)
%     tformAffine = pcregistercpd(moving,fixed,'Transform','Affine','Verbose',true);
%     ptCloudAffine = pctransform(moving, tformAffine);
% 
%     %% nonrigid registration (CPD) from original moving
%     tformNonrigid = pcregistercpd(moving,fixed,'Verbose',true);
%     ptCloudNonrigid = pctransform(moving, tformNonrigid);
%     
%     %% nonrigid registration (CPD)from non-rigid
%     tformNonrigid = pcregistercpd(ptCloudRigid,fixed,'Transform','nonrigid','Verbose',true);
%     ptCloudNonrigid = pctransform(ptCloudRigid, tformNonrigid);
    
end

    %% Visualization
%     axs = [];
%     figure; tiledlayout(2,3);
%     
%     % Before registration
%     nexttile;
%     pcshow(ptCloudRef); axis off;
%     axs = [axs, gca];
%     nexttile;
%     pcshow(ptCloudCurrent); axis off;
%     axs = [axs, gca];
%     nexttile;
%     pcshowpair(ptCloudRef,ptCloudCurrent); axis off;
%     axs = [axs, gca];
%     
%     % after rigid registration
%     nexttile;
%     pcshow(ptCloudRef); axis off;
%     axs = [axs, gca];
%     nexttile;
%     pcshow(ptCloudRigid); axis off;
%     axs = [axs, gca];
%     nexttile;
%     pcshowpair(ptCloudRef,ptCloudRigid); axis off;
%     axs = [axs, gca];
%     
%     hlink = linkprop(axs,{'CameraPosition','CameraUpVector','CameraViewAngle',...
%         'CameraViewAngleMode','CameraUpVectorMode','CLim','Colormap'}); % ,'View'

end

