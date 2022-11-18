function [distance,vnormal] = surfaceDisplace2Distance(vertex,triedge,displacement)
    %% convert surface displacement to distance
    
    %% Compute the vertex norm
    vnormal = vertexNormal(triangulation(double(triedge),double(vertex)));
    
    %% convert distance
    distance = transpose(dot(displacement,vnormal,2));
end
    