function [vertexReordered,pcLocSub] = vertexReorder(vertexTarget, vertexMean)
        %% reorder the points based on the closest points from the group mean (KNN)
        [pcLocSub,~] = knnsearch(vertexTarget, vertexMean);

        %% Recreate vertexReordered using the rearranged locations
        vertexReordered = vertexTarget(pcLocSub,:);
