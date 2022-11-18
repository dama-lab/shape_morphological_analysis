function [surfDistance, surfDispField] = vol2surfDisplacement(vertex, triedge, volDisplacement)
    %% find the surface displacement from a volume displacement field
    % Input
    %   vertexSub: vertex subscription index
    %   triedge:   triangle-edge of the surface connecting all the vertice
    %   volDisplacement: 4D displacement field
    
    %%
    surfDispField = [];
    fprintf(' interp ...');
    for dim = 1:size(volDisplacement,4)
        %% Loop through x/y/z
        fprintf(' dim %d ',dim);
        displacementDim = volDisplacement(:,:,:,dim);
        surfDispField(:,dim) = vol2vertexData(vertex,displacementDim,[],[],8);
    end
    
    %% Calculate surface displacement from displacement field
    surfDistance = surfaceDisplace2Distance(vertex,triedge,surfDispField);
    surfDistance = surfDistance';
end