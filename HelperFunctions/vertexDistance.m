function [distance,vertexReordered,displacement]  = vertexDistance(vertexTarget, vertexMean, vnormal)
    %% reorder the vertex
    vertexReordered = vertexReorder(vertexTarget, vertexMean);
    %% displacement (along x,y,z axis)
    displacement = vertexMean - vertexReordered;
    %% distance (scalar value, signed)
    distance = transpose(dot(displacement,vnormal,2));
end