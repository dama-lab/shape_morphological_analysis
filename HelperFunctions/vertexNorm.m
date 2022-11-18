function [vnormal, TR] = vertexNorm(triedge,vertex)
    %% Find the norm of displacement for each individual to the group mean
    % vnorm Ref: mial-tools/matlab/VolSubCortShapeAnalysis/generate_norm_disp (Ref to subcortical group difference module)
    
    TR = triangulation(double(triedge),double(vertex));
    vnormal = vertexNormal(TR);
end