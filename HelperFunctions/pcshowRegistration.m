function [fig,axs] = pcshowRegistration(ptCloudRef, ptCloudFloat, ptCloudRigid)
    %% Visualization effect of point cloud registratoin
    % ptCloudRef: Reference point cloud
    % ptCloudCurrent: original unregistered moving/floating point cloud
    % ptCloudRigid: registered point cloud
    
    axs = [];
    fig = figure;
    tiledlayout(2,3);
    
    % Before registration
    nexttile;
    pcshow(ptCloudRef,'MarkerSize',50); axis off;
    axs = [axs, gca];
    nexttile;
    pcshow(ptCloudFloat,'MarkerSize',50); axis off;
    axs = [axs, gca];
    nexttile;
    pcshowpair(ptCloudRef,ptCloudFloat); axis off;
    axs = [axs, gca];
    
    % after registration
    nexttile;
    pcshow(ptCloudRef,'MarkerSize',50); axis off;
    axs = [axs, gca];
    nexttile;
    pcshow(ptCloudRigid,'MarkerSize',50); axis off;
    axs = [axs, gca];
    nexttile;
    pcshowpair(ptCloudRef,ptCloudRigid); axis off;
    axs = [axs, gca];
%     
%     hlink = linkprop(axs,{'CameraPosition','CameraUpVector','CameraViewAngle',...
%         'CameraViewAngleMode','CameraUpVectorMode','Colormap'}); % ,'CLim' ,'View'
end
