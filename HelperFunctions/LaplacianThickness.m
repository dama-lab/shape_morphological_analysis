function [thickness, LalpacianField, L0, L1] = LaplacianThickness(volSeg,labelCentre,labelBoundary)
    %% Step 1: Calculate Laplacian field
    fprintf('\nLaplacian Field: ');
    LalpacianField = LaplaceField(volSeg, labelCentre, labelBoundary, 1);
    %% Step 2: Calculate Laplacian thickness
    fprintf('\nLaplacian Thickness: ');
    [thickness, L0, L1] = laplace_thickness(LalpacianField);