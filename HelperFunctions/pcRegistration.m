function [pcRegistered, tformRigid] = pcRegistration(pcFloat,pcRef)
%PCREGISTRATION 
%   Input:
%       pcMoving: 
%       pcFix
    %% downsample to speed up and improve the accuracy of the registration
    gridSize = 0.1;
    fixed  = pcdownsample(pcRef,'gridAverage', gridSize);
    moving = pcdownsample(pcFloat,'gridAverage', gridSize);

    %% rigid registration (ICP)
    [tformRigid, rmse] = pcregistericp(moving, fixed, 'Metric', 'pointToPoint', 'Extrapolate', true);
    pcRegistered = pctransform(pcFloat,tformRigid);

end

