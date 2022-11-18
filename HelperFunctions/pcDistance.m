function [distance, vertexTarget, displacement, pcTarget, pcLocSub] = pcDistance(pcTarget, pcMean, vnormal)

        %% Test point location using color code to determine point locations
        % pcTarget = pointCloud(pcTarget.Location,'Intensity',pcIdIntensity);
        % pcshow(pcTarget);
        % pcshowpair(pcMean, pcTarget);

        %% reorder the points based on the closest points from the group mean (KNN)
        [pcTarget, pcLocSub] = pcReorder(pcTarget, pcMean);
        
%         [pcLocSub,~] = knnsearch(pcTarget.Location,pcMean.Location);
% 
%         %% Recreate pcTarget using the rearranged locations
%         pcTarget = pointCloud(pcTarget.Location(pcLocSub,:));
%         %pcshow(pcTarget);
%         %pcTarget.Location = pcTarget.Location(pcLocSub,:);            


        %% For calculating groupwised average shape
        vertexTarget = pcTarget.Location;
        %% displacement (along x,y,z axis)
        displacement = (pcMean.Location - vertexTarget);
        %% distance (scalar value, signed)
        distance = transpose(dot(displacement,vnormal,2)); 
