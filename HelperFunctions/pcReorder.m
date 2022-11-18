function [pcTarget,pcLocSub] = pcReorder(pcTarget, pcMean)

    %% reorder the points based on the closest points from the group mean (KNN)
    [pcLocSub,~] = knnsearch(pcTarget.Location,pcMean.Location);

    %% Recreate pcTarget using the rearranged locations
    pcTarget = pointCloud(pcTarget.Location(pcLocSub,:));
    
    %pcshow(pcTarget);
    %pcTarget.Location = pcTarget.Location(pcLocSub,:);            

