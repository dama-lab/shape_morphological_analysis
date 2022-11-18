function [thickness, LalpacianField, L0, L1] = LaplacianThickness(volSeg,labelCentre,labelBoundary)
    %% labelCentre: label number for the structure of interest for calculating the thickness (e.g. gray matter)
    %% labelBoundary: label number of the adjacent structure with shared boundary of the structure of interest (e.g. white matter)
    
    %% Step 1: Calculate Laplacian field
    fprintf('\nLaplacian Field: ');
    LalpacianField = LaplaceField(volSeg, labelCentre, labelBoundary, 1);
    %% Step 2: Calculate Laplacian thickness
    fprintf('\nLaplacian Thickness: ');
    [thickness, L0, L1] = laplace_thickness(LalpacianField);