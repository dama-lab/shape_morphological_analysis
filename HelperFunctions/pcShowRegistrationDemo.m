

ptCloudRef = pointCloud(surfGroup('groupMean'));
ptCloudFloat = pcTarget;
ptCloudRigid = pcGroup(target);
[fig,axs] = pcshowRegistration(ptCloudRef, ptCloudFloat, ptCloudRigid);
hlink = linkprop(axs,{'CameraPosition','CameraUpVector','CameraViewAngle',...
        'CameraViewAngleMode','CameraUpVectorMode','Colormap'}); %,'CLim' % ,'View'